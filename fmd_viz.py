import os
import yaml
import json
import pandas as pd
import plotly.graph_objects as go
from pathlib import Path
from lxml import etree
import pint
from dash import Dash, dcc, html, Input, Output
import click
import logging

# --- CONFIGURATION ---
RESULTS_FOLDER = Path("./data/results")
u = pint.UnitRegistry()
NS = {'ipc': 'http://webstds.ipc.org/175x/2.0'}

def get_node_mass(node, ns):
    unit = node.get('UOM')
    value_str = node.get('value')
    if value_str is None:
        amount_node = node.find('ipc:Amount', ns)
        if amount_node is not None:
            value_str = amount_node.get('value')
            unit = amount_node.get('UOM')
        else:
            return 0 * u.gram
    try:
        return u.parse_expression(f"{value_str} {unit}").to('gram')
    except:
        return 0 * u.gram

def build_xml_links(node, parent_name, links):
    """
    Recursively traverses the XML to build parent->child links.
    Captures SubProduct, HomogeneousMaterial, and Substance.
    """
    # Look for children of types: SubProduct, HomogeneousMaterial, Substance
    targets = node.xpath('ipc:SubProduct | .//ipc:HomogeneousMaterial | .//ipc:Substance', namespaces=NS)
    
    for child in targets:
        # Determine Name
        if child.tag.endswith('SubProduct'):
            id_node = child.find(".//ipc:ProductID", NS)
            child_name = id_node.get("itemName") if id_node is not None else "Unknown SubProduct"
            child_type = 'SubProduct'
        else:
            child_name = child.get("name")
            child_type = 'HomogeneousMaterial' if child.tag.endswith('HomogeneousMaterial') else 'Substance'
        
        mass = get_node_mass(child if not child.tag.endswith('SubProduct') else child.find(".//ipc:ProductID", NS), NS)
        
        if mass > 0:
            links.append({
                'source': parent_name,
                'target': child_name,
                'value': mass,
                'type': child_type
            })
            # Recursive call for nested structures (only SubProducts usually have children)
            build_xml_links(child, child_name, links)

def treat_node_for_viz(node, parent_name, ipcs, ns, previous_node_type):
    """
    Builds a list of dicts for the Sankey DataFrame.
    parent_name: The string name of the parent node
    ipcs: List of XPaths for [SubProduct, HomogeneousMaterial, Substance]
    """
    links = []

    # Find children for the current search pattern (ipcs[0])
    children = node.findall(ipcs[0], ns)
    node_type = ipcs[0].split(':')[-1]

    for child in children:
        # 1. Safely extract the name as a STRING
        if node_type in ["SubProduct"]:
            id_node = child.find(".//ipc:ProductID", ns)
            # Use only the attribute name string, don't pass the ns dict to .get()
            name = id_node.get("itemName") if id_node is not None else "Unknown Product"
        else:
            name = child.get("name")
        
        # Fallback if name is missing or unexpectedly not a string
        name = str(name) if name is not None else "Unnamed"

        # 2. Extract Mass (reusing your get_mass logic)
        # Assuming get_mass is defined globally and returns a pint object
        try:
            mass_obj = get_node_mass(child if node_type not in ["SubProduct"] else child.find(".//ipc:ProductID", ns), ns)
            mass_value = mass_obj.to('gram').magnitude
        except:
            mass_value = 0.0

        # 3. Create the Link
        links.append({
            'source': str(parent_name),
            'target': name,
            'value': float(mass_value),
            'type': node_type

        })

        # 4. Recursion: If there are more patterns in ipcs, go deeper
        if len(ipcs) > 1:
            # Pass the child as the new parent for the next level
            links.extend(treat_node_for_viz(child, name, ipcs[1:], ns, previous_node_type=node_type))

    return links

def load_mapping_dfs(mappings_folder):
    """Load all mapping CSVs and normalize keys for fast lookup."""
    mapping_files = {
        'SubProduct': 'subproduct_mapping.csv',
        'HomogeneousMaterial': 'homo_material_mapping.csv',
        'Substance': 'substance_mapping.csv'
    }
    
    mapping_dfs = {}
    for node_type, filename in mapping_files.items():
        filepath = Path(mappings_folder) / filename
        try:
            df = pd.read_csv(filepath)
            df['normalized_key'] = df.iloc[:, 0].str.strip().str.lower()
            mapping_dfs[node_type] = df
        except FileNotFoundError:
            logging.warning(f"Mapping file not found: {filepath}")
            mapping_dfs[node_type] = pd.DataFrame()
    
    return mapping_dfs

def is_mapped(target_name, node_type, mapping_dfs):
    """Check if a target is mapped (has non-empty ecoinvent column)."""
    if node_type not in mapping_dfs or mapping_dfs[node_type].empty:
        return False
    
    normalized_target = str(target_name).strip().lower()
    df = mapping_dfs[node_type]
    
    # Find the row matching this target
    match = df[df['normalized_key'] == normalized_target]
    
    if match.empty:
        return False
    
    # Check if 'ecoinvent activity' column exists and is not empty
    if 'ecoinvent activity' in match.columns:
        return pd.notna(match.iloc[0]['ecoinvent activity']) and \
               str(match.iloc[0]['ecoinvent activity']).strip() != ''
    
    return False

def create_app(fmd_folder, mappings_folder):
    app = Dash(__name__)

    app.layout = html.Div([
        html.H2("FMD Sankey: Raw XML vs YAML Result", style={'fontFamily': 'sans-serif'}),
        
        html.Div([
            html.Div([
                html.Label("Source XML File:"),
                dcc.Dropdown(id='xml-selector', placeholder="Select XML...")
            ], style={'width': '30%', 'display': 'inline-block', 'marginRight': '20px'}),

            html.Div([
                html.Label("Mass Cutoff (percentage of product mass):"),
                dcc.Input(id='cutoff-input', type='number', value=1, step=1)
            ], style={'display': 'inline-block', 'marginRight': '20px'}),

            html.Div([
                html.Label("Color Mode:"),
                dcc.Dropdown(
                    id='color-mode',
                    options=[
                        {'label': 'Architecture', 'value': 'arch'},
                        {'label': 'Highlight Dictionary', 'value': 'dict'},
                        {'label': 'Mapping availability', 'value': 'mapping'}
                    ], value='arch'
                )
            ], style={'width': '200px', 'display': 'inline-block'}),
            
            # New Textarea for the Highlight Dictionary
            html.Div(id='dict-input-container', children=[
                html.Label("Highlight Dictionary (JSON format):"),
                dcc.Textarea(
                    id='highlight-dict',
                    value='{"Copper":"red"}',
                    placeholder='{"Gold": "gold", "Copper": "#b87333"}',
                    style={'width': '100%', 'height': '60px', 'fontFamily': 'monospace'}
                )
            ], style={'marginTop': '10px', 'display': 'none'})
            
        ], style={'padding': '20px', 'background': '#f9f9f9', 'borderRadius': '8px','flex': '0 1 auto'}),
        
        html.Div(
                children=[
                    dcc.Graph(
                        id='sankey-graph',
                        responsive=False,
                        config={'responsive': False, 'displayModeBar': True},
                        style={'width': '100%', 'height': '100%'}
                    ),
                ], 
                style={
                    'flex': '1 1 auto',
                    'overflowY': 'auto',
                    'overflowX': 'hidden',
                    'border': '1px solid #ccc',
                    'position': 'relative',
                    'minHeight': '0',
                    'width': '100%'
                }
            ),
        dcc.Interval(id='refresh', interval=5000)
    ], style={
        'display': 'flex', 
        'flexDirection': 'column', 
        'height': '100%',
        'width': '100%',
        'margin': '0', 
        'padding': '0',
        'overflow': 'hidden'
    })

    @app.callback(
        Output('dict-input-container', 'style'),
        Input('color-mode', 'value')
    )
    def toggle_dict_input(mode):
        if mode == 'dict':
            return {'marginTop': '10px', 'display': 'block'}
        return {'marginTop': '10px', 'display': 'none'}

    @app.callback(
        Output('xml-selector', 'options'),
        [Input('refresh', 'n_intervals')]
    )
    def update_file_dropdown(n):
        fmd_path = Path(fmd_folder)
        xml_files = []
        
        if fmd_path.exists():
            xml_files.extend(list(fmd_path.glob("*.xml")))
        
        base_path = Path(__file__).parent
        xml_files.extend(list(base_path.glob("*.xml")))
        
        return [{'label': f.name, 'value': str(f.absolute())} for f in set(xml_files)]
    @app.callback(
        Output('cutoff-input', 'value'),
        Input('color-mode', 'value'),
        prevent_initial_call=True
    )
    def reset_cutoff_on_mode_change(mode):
        return 1
    @app.callback(
        Output('sankey-graph', 'figure'),
        [Input('xml-selector', 'value'),
        Input('cutoff-input', 'value'),
        Input('color-mode', 'value'),
        Input('highlight-dict', 'value')]
    )

    def update_graph(xml_path, cutoff, mode, dict_str):
        if not xml_path or not Path(xml_path).exists():
            return go.Figure().update_layout(title="Select an XML file to visualize raw architecture.")

        # Parse XML
        parser = etree.XMLParser(ns_clean=True)
        tree = etree.parse(xml_path, parser)
        
        product_node = tree.find('.//ipc:Product', NS)
        prod_id = product_node.find('ipc:ProductID', NS)
        root_name = prod_id.get('itemName')
        product_mass_grams = get_node_mass(prod_id, NS).to('gram').magnitude
        
        ipcs_patterns = ['ipc:SubProduct', './/ipc:HomogeneousMaterial', './/ipc:Substance']    
        links_data = treat_node_for_viz(product_node, root_name, ipcs_patterns, NS, "Product")
        
        full_df = pd.DataFrame(links_data)
        if full_df.empty:
            return go.Figure().update_layout(title="No data found in XML file.")

        full_df['value'] = pd.to_numeric(full_df['value'], errors='coerce')
        cutoff_val = cutoff if cutoff is not None else 1
        threshold = float(product_mass_grams * (cutoff_val / 100.0))
        full_df = full_df[full_df['value'] >= threshold]

        # --- STEP 1: SPLIT AND AGGREGATE BY TYPE ---
        df_sp = full_df[full_df['type'] == 'SubProduct'].groupby(['source', 'target', 'type'], as_index=False)['value'].sum()
        df_hm = full_df[full_df['type'] == 'HomogeneousMaterial'].groupby(['source', 'target', 'type'], as_index=False)['value'].sum()
        df_sub = full_df[full_df['type'] == 'Substance'].groupby(['source', 'target', 'type'], as_index=False)['value'].sum()

        # Combine back into a clean dataframe
        df = pd.concat([df_sp, df_hm, df_sub], ignore_index=True)

        # --- STEP 2: CREATE UNIQUE NODE IDENTIFIERS ---
        def get_node_id(name, n_type):
            if name == root_name: return f"{name}_Product"
            return f"{name}_{n_type}"

        source_type_map = {row['target']: row['type'] for _, row in df.iterrows()}
        source_type_map[root_name] = 'Product'

        def get_source_type(row):
            if row['source'] == root_name: return 'Product'
            if row['type'] == 'Substance': return 'HomogeneousMaterial'
            if row['type'] == 'HomogeneousMaterial': return 'SubProduct'
            return 'SubProduct'

        df['source_id'] = df.apply(lambda r: get_node_id(r['source'], get_source_type(r)), axis=1)
        df['target_id'] = df.apply(lambda r: get_node_id(r['target'], r['type']), axis=1)

        # Build the final node list for the Sankey
        unique_nodes = pd.concat([
            df[['source', 'source_id']].rename(columns={'source': 'label', 'source_id': 'id'}),
            df[['target', 'target_id']].rename(columns={'target': 'label', 'target_id': 'id'})
        ]).drop_duplicates('id')

        nodes_labels = unique_nodes['label'].tolist()
        nodes_ids = unique_nodes['id'].tolist()
        n_idx = {nid: i for i, nid in enumerate(nodes_ids)}

        # --- STEP 3: ASSIGN COLUMNS (X COORDS) ---
        type_to_x = {'Product': 0.01, 'SubProduct': 0.33, 'HomogeneousMaterial': 0.66, 'Substance': 0.99}
        node_x = []
        for nid in nodes_ids:
            ntype = nid.split('_')[-1]
            node_x.append(type_to_x.get(ntype, 0.5))

        # --- STEP 4: COLORING ---
        architecture_c_map = {'Product': '#3b2c7d', 'SubProduct': '#446e9b', 'HomogeneousMaterial': '#d99a29', 'Substance': '#3cb371'}
        bg_node = "#E5E5E5"

        if mode == 'mapping':
            mapping_dfs = load_mapping_dfs(mappings_folder)
            node_colors = []
            for _, row in unique_nodes.iterrows():
                ntype = row['id'].split('_')[-1]
                if is_mapped(row['label'], ntype, mapping_dfs):
                    node_colors.append(architecture_c_map.get(ntype, '#999'))
                else:
                    node_colors.append(bg_node)
            link_colors = [architecture_c_map.get(row['type'], '#ccc') if is_mapped(row['target'], row['type'], mapping_dfs) else "rgba(200, 200, 200, 0.2)" for _, row in df.iterrows()]
        
        elif mode == 'dict':
            user_map = {}
            if dict_str:
                try: user_map = json.loads(dict_str)
                except: pass
            node_colors = [next((color for key, color in user_map.items() if key in name), bg_node) for name in nodes_labels]
            link_colors = [user_map.get(row['target'], "#F0F0F0") for _, row in df.iterrows()]
        
        else: # Architecture Mode
            node_colors = [architecture_c_map.get(nid.split('_')[-1], '#999') for nid in nodes_ids]
            link_colors = [architecture_c_map.get(t, '#ccc') for t in df['type']]

        # --- STEP 5: BUILD FIGURE ---
        fig = go.Figure(go.Sankey(
            arrangement='snap',
            node=dict(
                pad=15, thickness=20,
                label=nodes_labels,
                color=node_colors,
                x=node_x,
                y=[0.1] * len(nodes_ids),
                hovertemplate='Node: %{label}<br>Mass: %{value:.4f}g<extra></extra>'
            ),
            link=dict(
                source=df['source_id'].map(n_idx),
                target=df['target_id'].map(n_idx),
                value=df['value'],
                color=link_colors,
                hovertemplate='Path: %{source.label} → %{target.label}<br>Mass: %{value:.4f}g<extra></extra>'
            )
        ))

        if mode == 'arch':
            legend_data = [(key,color) for key, color in architecture_c_map.items()]
            legend_title = "Architecture Legend"
        elif mode == 'mapping':
            legend_data = [(key,color) for key, color in architecture_c_map.items()]
            legend_data.append(("Unmapped node", bg_node))
            legend_title = "Architecture Legend"
        else:
            user_map = {}
            if dict_str:
                try:
                    user_map = json.loads(dict_str)
                except:
                    pass
            
            legend_data = [(f"{k}", v) for k, v in user_map.items()]
            legend_data.append(("Others", bg_node))
            legend_title = "Highlight Legend"
        
        for label, color in legend_data:
            fig.add_trace(go.Scatter(
                x=[None], y=[None],
                mode='markers',
                marker=dict(size=12, color=color, symbol='square'),
                showlegend=True,
                name=label
            ))

        # --- CALCULATE DYNAMIC HEIGHT BY NODE TYPE ---
        type_heights = {}
        threshold_percent = 2.0  # 2% threshold
        threshold_mass = product_mass_grams * (threshold_percent / 100.0)
        product_height = 50  # Fixed height for root node
        px_per_percent = 0.1

        for node_type in ['SubProduct', 'HomogeneousMaterial', 'Substance']:
            type_nodes = unique_nodes[unique_nodes['id'].str.endswith(node_type)]
            height = 0
            
            for _, node in type_nodes.iterrows():
                # Find all links targeting this node
                node_mass = df[df['target_id'] == node['id']]['value'].sum()
                mass_percent = (node_mass / product_mass_grams) * 100

                if node_mass >= threshold_mass:
                    # For nodes > 2%: allocate height proportional to mass contribution
                    height += mass_percent * px_per_percent
                else:
                    # For nodes <= 2%: fixed 30px each
                    height += 30
            
            type_heights[node_type] = height
        
        max_height = max(max(type_heights.values()), product_height) + 200  # Add buffer for margins and legend
        dynamic_height = max(600, max_height)

        fig.update_layout(
            title_text=f"Visualizing: {Path(xml_path).name}",
            height=dynamic_height,
            margin=dict(t=50, b=20, l=10, r=10),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, visible=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, visible=False),
            legend=dict(
                title=legend_title,
                orientation="h",
                yanchor="bottom",
                y=1.02,
                xanchor="right",
                x=1
            ),
        )
        
        return fig
    return app
@click.command()
@click.option("-f", "--fmd_folder", default="./FMD", help="FMD folder for XML files")
@click.option("-m", "--mappings_folder", default="./data", help="Folder containing the mappings")
@click.option("-v", "--verbose", count=True, help="Increase verbosity (-v, -vv, -vvv)")
def run_viz(fmd_folder, mappings_folder, verbose):
    """
    Run the FMD Sankey visualization web app.
    """
    
    level = logging.WARNING
    if verbose == 1:
        level = logging.INFO
    elif verbose >= 2:
        level = logging.DEBUG
    
    logging.basicConfig(level=level)
    
    logging.info(f"Using FMD folder: {fmd_folder}")
    logging.info(f"Using mappings folder: {mappings_folder}")
    
    app = create_app(fmd_folder, mappings_folder)
    app.run(debug=True, port=8050)

if __name__ == "__main__":
    run_viz()