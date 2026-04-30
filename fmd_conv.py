import click
from lxml import etree
import numpy as np
from pathlib import Path
import pandas as pd
import re
import yaml

def load_mapping(map_folder):
    def df_to_map(path, key_col):
        df = pd.read_csv(path)
        df[key_col] = df[key_col].str.strip().str.lower()
        return df.set_index(key_col)[["ecoinvent activity", "location"]] \
                 .apply(tuple, axis=1) \
                 .to_dict()

    mapping_dict = {
        **df_to_map(f"{map_folder}/substance_mapping.csv", "substance name"),
        **df_to_map(f"{map_folder}/homo_material_mapping.csv", "homogeneous material"),
    }

    return lambda name: mapping_dict[name.strip().lower()]

default_ac = {
    "Substance": ('market for copper, anode', 'GLO'),
    "HM": ('metal working, average for copper product manufacturing', 'RoW')
    }

def treat_node(node, sp_name, hm_name, ns, apply_mapping, node_type):
    this_name = node.get('name', 'Unknown')

    act_name, location = apply_mapping(this_name)

    value_str = node.get('value') # Try to get value directly from the Substance node
    if value_str is None: # If not present, look for a child <Amount> node
        amount_node = node.find('ipc:Amount', ns)
        if amount_node is not None:
            value_str = amount_node.get('value', '0')
        else:
            value_str = '0'
            print(f"Error: The value for the homogeneous {this_name} in file {xml_file.name} is zero")
    value = float(value_str)

    if pd.isna(act_name):
        act_name, location = default_ac[node_type]
        value = 0

    return {
        'act_name': act_name,
        'location': location,
        'c_subproduct': sp_name,
        'c_homogeneous_material': hm_name,
        'c_substance': 'process' if type == "HM" else this_name,
        'amount': {
            'value': value,
            'unit': 'mg'
        }
    }

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
    
    # Load XML
    parser = etree.XMLParser(ns_clean=True)
    tree = etree.parse(xml_file, parser)
    
    # Define namespaces
    ns = {'ipc': 'http://webstds.ipc.org/175x/2.0'}

    # Get product info
    product_node = tree.find('.//ipc:Product', ns)
    product_id_node = product_node.find('ipc:ProductID', ns)
    product_name = f"{product_id_node.get('itemName')}_{product_id_node.get('version')}"
    total_mass = float((product_id_node.find('ipc:Amount', ns).get('value', 0)) if product_id_node.find('ipc:Amount', ns) is not None else 0)
    
    inputs = {}
    computed_mass = 0.0

    # ------------------------------------------------------
    subproducts = product_node.findall('ipc:SubProduct', ns)
    for sp in subproducts: # iteration in all subproducts
        sp_product_node = sp.find('ipc:ProductID', ns)
        sp_name = sp_product_node.get('itemName', f"{len(inputs):03d}") \
            if sp_product_node is not None else f"{len(inputs):03d}"

        hmlist = sp.findall('.//ipc:HomogeneousMaterial', ns)
        for hm in hmlist:
            hm_name = hm.get('name', f"HM_{len(inputs)}")
            inputs[f"process_input_{len(inputs):03d}"] = treat_node(hm, sp_name, hm_name, ns, apply_mapping, 'HM')

            subs = hm.findall('.//ipc:Substance', ns) 
            
            for sub in subs:
                new_sub = treat_node(sub, sp_name, hm_name, ns, apply_mapping, 'Substance')
                inputs[f"process_input_{len(inputs):03d}"] = new_sub
                computed_mass += new_sub["amount"]["value"]

    tolerance = 0.01  # mg
    if abs(computed_mass - total_mass) > tolerance:
        print(f"WARNING {xml_file}: mismatch total_mass={total_mass:.2f} vs sum_substances={computed_mass:.2f}")

    yaml_data = {
        'output': {
            'product': product_name,
            'amount': {'value': 1, 'unit': 'unit'}
        },
        'c_total_mass': total_mass,
        'c_accounted_mass': computed_mass,
    }

    yaml_data['inputs'] = inputs

    yaml_file = output_folder / f"{product_name}.yaml"
    with open(yaml_file, 'w') as f:
        yaml.dump(yaml_data, f, sort_keys=False)


@click.command()
@click.argument("input_files", nargs=-1, type=click.Path(exists=True))
@click.option("-o", "--output_folder", default="./data/results", help="Output folder for results")
@click.option("-m", "--mappings_folder", default="./data", help="Folder containing the mappings")
def run_conv(input_files, output_folder, mappings_folder):
    """
    Translate to YAML one or multiple XML FMD files.
    """

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
