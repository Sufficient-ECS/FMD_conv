import click
from lxml import etree
from pathlib import Path
import pandas as pd
import re
import yaml
import logging

def update_mappings(files, output_folder):
    ns = {'ipc': 'http://webstds.ipc.org/175x/2.0'}

    file_names = ["subproduct", "homo_material", "substance"]
    c_names = ["subproduct name", "homogeneous material", "substance name"]
    ipcs = [".//ipc:SubProduct/ipc:ProductID", ".//ipc:HomogeneousMaterial", ".//ipc:Substance"]
    name_cells = ["itemName", "name", "name"]

    ds = [pd.read_csv(output_folder / f"{file_name}_mapping.csv")\
            .assign(**{c_name: lambda df: df[c_name].str.lower()}) \
            .set_index(c_name)\
            .to_dict(orient="index")
            for file_name, c_name in zip(file_names, c_names)]

    for file in files:
        try:
            tree = etree.parse(file)
            root = tree.getroot()

            for i, d in enumerate(ds):
                for sub in root.findall(ipcs[i], ns):
                    name = sub.get(name_cells[i]).strip().lower()
                    if name not in d:
                        logging.info(f"Adding {file_names[i]} {name}")
                        d[name] = {"ecoinvent activity": None, "location": None}

        except Exception as e:
            print(f"Error processing {file}: {e}")

    for file_name, c_name, d in zip(file_names, c_names, ds):
        pd.DataFrame.from_dict(d, orient="index") \
          .rename_axis(c_name) \
          .sort_index() \
          .to_csv(output_folder / f"{file_name}_mapping.csv")

@click.command()
@click.argument("input_files", nargs=-1, type=click.Path(exists=True))
@click.option("-o", "--output_folder", default="./data", help="Output folder for results")
@click.option("-v", "--verbose", count=True, help="Increase verbosity (-v, -vv, -vvv)")
def run_conv(input_files, output_folder, verbose):
    """
    Run LCA impacts on one or multiple YAML foreground files.
    """

    level = logging.WARNING  # default
    if verbose == 1:
        level = logging.INFO
    elif verbose >= 2:
        level = logging.DEBUG

    logging.basicConfig(level=level)

    if not input_files:
        raise click.UsageError("You must provide at least one input file.")
    
    # Ensure output directory exists
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    update_mappings(input_files, output_folder)


if __name__ == "__main__":
    run_conv()
