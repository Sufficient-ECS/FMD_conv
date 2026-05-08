import click
import logging
from lxml import etree
import numpy as np
from pathlib import Path
import pandas as pd
import pint
import re
import yaml

u = pint.UnitRegistry()

def load_mapping(map_folder):
    def df_to_map(path, key_col):
        df = pd.read_csv(path)
        df[key_col] = df[key_col].str.strip().str.lower()
        return df.set_index(key_col)[["ecoinvent activity", "location", "process"]] \
                 .apply(tuple, axis=1) \
                 .to_dict()

    mapping_dict = {
        **df_to_map(f"{map_folder}/substance_mapping.csv", "substance name"),
        **df_to_map(f"{map_folder}/homo_material_mapping.csv", "homogeneous material"),
        **df_to_map(f"{map_folder}/subproduct_mapping.csv", "subproduct name"),
    }
    
    return lambda name: mapping_dict[name.strip().lower()]

def get_mass(node, ns):
    unit = node.get('UOM')
    value_str = node.get('value') # Try to get value directly from the Substance node
    if value_str is None: # If not present, look for a child <Amount> node
        amount_node = node.find('ipc:Amount', ns)
        if amount_node is not None:
            value_str = amount_node.get('value')
            unit = amount_node.get('UOM')
        else:
            return 0 * u.gram
    return u.parse_expression(f"{value_str} {unit}")

def treat_node(node, inputs, names, ipcs, ns, apply_mapping, prev_mass = None):
    computed_mass = 0 * u.gram

    for i in node.findall(ipcs[0], ns):
        this_names = names.copy()
        ind = this_names.index("_")


        if ind == 0: # If is first layer
            id_node = i.find(".//ipc:ProductID", ns)
            name = id_node.get("itemName", ns)
            mass = get_mass(id_node, ns)
        else:
            name = i.get("name", ns)
            mass = get_mass(i, ns)

        this_names[ind] = name

        map_info = apply_mapping(name)
        act_name, location = map_info[:2]
        if len(map_info) == 3:
            process = map_info[2]

        indent = (1 + ind) * "\t"
        is_accounted = not pd.isna(act_name)
        indicator = '!' if not is_accounted else ''
        if prev_mass != None: # If is not the first layer
            logging.debug(f"{indicator}{indent}{name} {mass} ({(mass/prev_mass).to('%'):.2f})")
        else:
            logging.debug(f"{indicator}{indent}{name} {mass}")

        if is_accounted:
            in_name = f"process_input_{len(inputs):03d}"
            inputs[in_name] =  {
                'act_name': act_name,
                'location': location,
                'c_subproduct': this_names[0],
                'c_homogeneous_material': this_names[1],
                'c_substance': this_names[2],
                'amount': {
                    'value': mass.magnitude,
                    'unit': str(mass.units)
                }
            }

        if this_names[-1] == "_" and (process or pd.isna(process)):
            computed_mass += treat_node(i, inputs, this_names, ipcs[1:], ns, apply_mapping, prev_mass = mass)
        else:
            computed_mass += mass

    return computed_mass

def sanitize_filename(name: str) -> str:
    # Replace invalid characters with a space
    name = re.sub(r'[<>:"/\\|?*\n\r\t]', ' ', name)
    # Collapse multiple spaces into one
    name = re.sub(r'\s+', ' ', name)
    return name.strip()

def ipc1752_to_yaml(xml_file: str, output_folder: str, apply_mapping):
    """
    Parse a single IPC-1752 XML BOM file into a flattened YAML activity.
    
    Rules:
    - One YAML file per XML (product)
    - Flatten SubProducts → HomogeneousMaterials → Substances
    - Keep unit as mg
    - Assign productref_input_001, _002, etc.
    - + homogeneous materials
    """

    logging.debug(xml_file)

    # Load XML
    parser = etree.XMLParser(ns_clean=True)
    tree = etree.parse(xml_file, parser)

    # Define namespaces
    ns = {'ipc': 'http://webstds.ipc.org/175x/2.0'}

    # Get product info
    product_node = tree.find('.//ipc:Product', ns)
    product_id_node = product_node.find('ipc:ProductID', ns)
    product_name = f"{product_id_node.get('itemName')}_{product_id_node.get('version')}"
    total_mass = get_mass(product_id_node, ns)

    logging.debug(product_name)

    inputs = {}
    computed_mass = 0.0 * u.gram

    ipcs = ['ipc:SubProduct', './/ipc:HomogeneousMaterial', './/ipc:Substance']
    computed_mass = treat_node(product_node, inputs, ["_", "_", "_"], ipcs, ns, apply_mapping)

    tolerance = 0.01 * u.mg
    if abs(computed_mass - total_mass) > tolerance:
        logging.warning(f"{xml_file}: mismatch total_mass={total_mass.to('mg'):.2f} vs sum_substances={computed_mass.to('mg'):.2f}")

    yaml_data = {
        'output': {
            'product': product_name,
            'amount': {'value': 1, 'unit': 'unit'}
        },
        'c_total_mass': str(total_mass),
        'c_accounted_mass': str(computed_mass),
        'c_item_number': product_id_node.get('itemNumber'),
    }

    yaml_data['inputs'] = inputs

    yaml_file = output_folder / f"{sanitize_filename(product_name)}.yaml"
    with open(yaml_file, 'w') as f:
        yaml.dump(yaml_data, f, sort_keys=False)


@click.command()
@click.argument("input_files", nargs=-1, type=click.Path(exists=True))
@click.option("-o", "--output_folder", default="./data/results", help="Output folder for results")
@click.option("-m", "--mappings_folder", default="./data", help="Folder containing the mappings")
@click.option("-v", "--verbose", count=True, help="Increase verbosity (-v, -vv, -vvv)")
def run_conv(input_files, output_folder, mappings_folder, verbose):
    """
    Translate to YAML one or multiple XML FMD files.
    """

    level = logging.WARNING  # default
    if verbose == 1:
        level = logging.INFO
    elif verbose >= 2:
        level = logging.DEBUG

    logging.basicConfig(level=level)

    if not input_files:
        raise click.UsageError("You must provide at least one input file.")

    # Write YAML
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    apply_mapping = load_mapping(mappings_folder)
    for xml_file in input_files:
        ipc1752_to_yaml(xml_file, output_folder, apply_mapping)


if __name__ == "__main__":
    run_conv()
