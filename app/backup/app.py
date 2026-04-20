import dash
from dash import Dash, html, dcc, Input, Output, State, no_update
import dash_bootstrap_components as dbc
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import seaborn as sns
from dash import dash_table
import statsmodels.api as sm
import matplotlib.colors as mcolors
import sys
import os
from .about import layout as about_page_layout

# Get the root directory (two levels up from the app.py file)
root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# Add the root directory to sys.path
sys.path.append(root_dir)

from sqlitedb.db_handler import DBHandler

# =========================
# Constants
# =========================

ALPHA = 0.05
LOESS_FRAC = 0.3

PHEWAS_HEIGHT = 600
PREVALENCE_HEIGHT = 300

Y_CAP_PERCENTILE = 99

# =========================
# Initialize DB and App
# =========================

# Create a connection to the SQLite database
database_file = "db/polygenie.db"
db_handler = DBHandler(database_file)

# Create the Dash app
app = Dash(__name__, external_stylesheets=[dbc.themes.FLATLY, '/assets/style.css'], suppress_callback_exceptions=True)
server = app.server

app.title = "PolyGenie"
print("App title is:", app.title)

# =========================
# Palette and Controls
# =========================

# Color palette definition
palette = sns.color_palette("colorblind")
colorblind_palette_hex = [mcolors.to_hex(color) for color in palette]

# Color palette graphs (Male, Female, All)
colors = ["#90BAAD", "#FF6542", "#56667A"]

# Create dropdown options
gwas_names = db_handler.get_gwas_names()
disease_options = [{'label': d, 'value': d} for d in gwas_names]

reference_options = [
    {'label': 'Bottom', 'value': 'low'},
    {'label': 'Bottom and Middle', 'value': 'rest'},
]

# Default division options (will be updated based on PRS when possible)
division_options = [{'label': '4 (quartiles)', 'value': '4'},
                    {'label': '10 (deciles)', 'value': '10'}]

# =========================
# Targe classes for tabs
# =========================

# Build dynamic tabs from the DB's target classes
_tc_df = db_handler.get_target_classes()

target_classes = (
    _tc_df['target_class'].astype(str).tolist()
    if not _tc_df.empty
    else ['Phecodes','ICD_codes','Metabolites','Questionnaire']
)

tabs_children = [dcc.Tab(label=tc, value=tc) for tc in target_classes]
default_tab = target_classes[0] if target_classes else 'ICD_codes'


col_headers = ['GWAS_code', 'Code', 'Reference', 'Division', 'Beta', 'OR', 'CI_5', 
               'CI_95', 'P', 'R2', 'logpxdir', 'Description', 'Domain', 'Class', 'Type']
col_headers_tokeep = ['GWAS_code', 'Code', 'Description', 'Domain', 'Beta', 'OR', 'CI_5', 'CI_95', 
                      'P']

def add_loess(fig, df, x, y, label, visible):
    """Add a LOESS-smoothed line to a Plotly figure."""
    if df.empty:
        return

    loess = sm.nonparametric.lowess(df[y], df[x], frac=LOESS_FRAC)

    fig.add_trace(
        go.Scatter(
            x=loess[:, 0],
            y=loess[:, 1],
            mode='lines',
            name=f'LOESS {label}',
            line=dict(width=3),
            visible=visible,
            showlegend=True
        )
    )

#####################################################################################################################
##################################################### Callbacks #####################################################
#####################################################################################################################

# Callback to update the graph based on the selected tab and dropdowns
@app.callback(
    Output('correlations-graph', 'figure'),
    [Input('disease-dropdown', 'value'),
     Input('reference-dropdown', 'value'),
     Input('division-dropdown', 'value'),
     Input('tabs', 'value')],
)
def update_graph(disease_value, reference_value, division_value, tab):
    """
    Callback function to update the graph contents

    :param disease_value: condition selected on the disease field of the filters.
    :param quartile_value: quartile selected on the quartile field of the filters.
    :param reference_value: reference selected on the reference field of the filters.
    :param division_value: division selected on the division field of the filters.
    :param tab: inidcates the active tab (to know which graph has to be shown).
    :return: figure with the updated graph content
    """
    
    target_class = tab

    if reference_value == "low + intermediate": reference_value = "rest"

    gwas = db_handler.get_gwas_code_from_name(disease_value)

    filtered_data = db_handler.get_correlations(
        gwas,
        reference_value,
        division_value,
        target_class
    )
    filtered_data = filtered_data.rename(columns={'target':'code'})
    filtered_data = filtered_data.dropna(subset=['P'])
    filtered_data = filtered_data.sort_values(by=['class', 'logpxdir'], ascending=[True, True])

    # Plotly figure based on filtered data
    x = 'description'
    y = 'logpxdir'
    title = (
    '' # XFR: I removed the title since it does not provide any useful info
        #'Metabolites vs logpxdir' if tab == 'met' 
        #else 'Variable vs. logpxdir' if tab == 'quest' 
        #else 'Diseases vs. logpxdir'
    )
    color = 'domain'

    unique_categories = filtered_data['domain'].unique()
    color_map = {category: colorblind_palette_hex[i % len(colorblind_palette_hex)] for i, category in enumerate(unique_categories)}
 
    filtered_data['odds_ratio_display'] = filtered_data['odds_ratio'].apply(
        lambda x: f"{x:.2f}" if pd.notna(x) else "N/A"
    )
    filtered_data = filtered_data.rename(columns={'target':'code'})
    fig = px.scatter(filtered_data, x=x, y=y,
                title=title,
                color= color,
                color_discrete_map=color_map,
                template="plotly_white",
                custom_data=['Code'],
                hover_data={
                    'Code': True, # Show internal code in hover
                    'description': True,  # Show the description
                    'class': True,
                    'beta': ':.4f',
                    'odds_ratio_display': True,
                    'odds_ratio': False,
                    'logpxdir': False,
                    'P': ':.2e'
                }
            )
    
    # Get the number of rows
    element_count = filtered_data['domain'].nunique()

    if element_count:

        positive_sig = -np.log10(0.05/element_count)
        positive_sig_text = f"{positive_sig:.3f}"

        negative_sig = np.log10(0.05/element_count)
        negative_sig_text = f"{negative_sig:.3f}"

        # Add the positive significance line
        fig.add_shape(
            type="line",
            x0=0,
            x1=1,
            xref='paper',
            y0=positive_sig,
            y1=positive_sig,
            yref='y',
            line=dict(
                color="red",
                width=2,
                dash="dash",
            ),
        )

        # Add a label to the significance line
        fig.add_annotation(
            xref="paper", 
            x=1, 
            y=positive_sig, 
            text=positive_sig_text, 
            showarrow=False, 
            font=dict(
                color="red",
                size=12
            ),
            align="right",
            yshift=10
        )

        # Add the negative significance line
        fig.add_shape(
            type="line",
            x0=0,
            x1=1,
            xref='paper',
            y0=negative_sig,
            y1=negative_sig,
            yref='y',
            line=dict(
                color="red",
                width=2,
                dash="dash",
            ),
        )

        # Add a label to the significance line
        fig.add_annotation(
            xref="paper", 
            x=1, 
            y=negative_sig, 
            text=negative_sig_text, 
            showarrow=False, 
            font=dict(
                color="red",
                size=12
            ),
            align="right",
            yshift=10
        )

        
        # Filter data points outside thresholds
        outliers = filtered_data[(filtered_data['logpxdir'] < negative_sig) | (filtered_data['logpxdir'] > positive_sig)]

        # Add annotations for outliers
        for index, row in outliers.iterrows():
            fig.add_annotation(
                x=row[x], 
                y=row['logpxdir'],
                text=f"{row[x]}",
                showarrow=False,
                font=dict(
                    color="black",
                    size=12
                ),
                align="center",
                yshift=10
            )
        
    xname = target_class.replace('_', ' ')

    # Update layout
    fig.update_layout(
        xaxis_title=xname,
        yaxis_title='log10(p-value) × Effect size direction',
        xaxis_visible=False,
        xaxis_showticklabels=False,
        showlegend=False
    ),

    # Add a black line at y=0
    fig.add_shape(
        type="line",
        x0=0,
        x1=1,
        xref='paper',
        y0=0,
        y1=0,
        yref='y',
        line=dict(
            color="black",
            width=1,
            ),
    )
    return fig

# Callback to update the prs table based on the selected filters
@app.callback(
    [Output('prs-table-content', 'data'),  # Dash DataTable expects 'data' as list of dicts
     Output('table-stored-data', 'data')],
    [Input('disease-dropdown', 'value'),
     Input('reference-dropdown', 'value'),
     Input('division-dropdown', 'value'),
     Input('tabs', 'value')]
)
def update_table(disease_value, reference_value, division_value, tab):
    """
    Callback function to update the content of the interactive DataTable and store filtered data.
    """

    target_class = tab

    if reference_value == "low + intermediate": reference_value = "rest"

    gwas = db_handler.get_gwas_code_from_name(disease_value)

    filtered_data = db_handler.get_correlations(
        gwas,
        reference_value,
        division_value,
        target_class
    )
    filtered_data = filtered_data.rename(columns={'target':'code'})

    filtered_data = filtered_data.dropna(subset=['P'])
    # Round numeric columns to 6 decimal places
    # Round numeric columns EXCEPT P
    numeric_columns = filtered_data.select_dtypes(include=[np.number]).columns
    numeric_columns = [c for c in numeric_columns if c != 'P']

    filtered_data[numeric_columns] = filtered_data[numeric_columns].round(6)

    # Ensure all data types are JSON serializable
    filtered_data = filtered_data.replace({np.nan: None})  # Replace NaN values with None for JSON compatibility
    # assign consistent column names used by the app
    filtered_data.columns = col_headers
    # select columns we want to show
    filtered_data = filtered_data[col_headers_tokeep]  # TODO: Need a more elegant solution

    # Sort by numeric P (if present) while values are still numeric, then format numeric columns for display
    if 'P' in filtered_data.columns:
        filtered_data = filtered_data.sort_values('P')

    # Format numeric display
    for c in ['Beta', 'OR', 'CI_5', 'CI_95']:
        if c in filtered_data.columns:
            filtered_data[c] = filtered_data[c].apply(lambda x: f"{x:.4f}" if pd.notnull(x) else None)

    # Format P-value into scientific notation with 2 decimal places for display
    if 'P' in filtered_data.columns:
        filtered_data['P'] = filtered_data['P'].apply(lambda x: f"{x:.2e}" if pd.notnull(x) else None)

    # Convert DataFrame to list of dictionaries for DataTable
    data_store = filtered_data.to_dict('records')

    return data_store, data_store  # Populate dash_table.DataTable and store data

# Callback to handle file download
@app.callback(
    Output("download-dataframe-excel", "data"),
    Input("download-button", "n_clicks"),
    State("table-stored-data", "data"),
    prevent_initial_call=True
)
def download_table(n_clicks, stored_data):
    """
    #Callback to download the table as an Excel file.
    """
    if not stored_data:
        return None

    # Convert stored data (list of dicts) back to DataFrame
    df = pd.DataFrame(stored_data)

    def to_excel(bytes_io):
        with pd.ExcelWriter(bytes_io, engine="xlsxwriter") as writer:
            df.to_excel(writer, sheet_name="Sheet1", index=False)
    return dcc.send_bytes(to_excel, "polygenie-table_data.xlsx")

def filter_values(unfiltered_data, disease_value, quartile_value, reference_value, division_value):
    """
    Function to filter the data according to the content of the filters

    :param unfiletered_data: original dataset
    :param disease_value: condition selected on the disease field of the filters.
    :param quartile_value: quartile selected on the quartile field of the filters.
    :param reference_value: reference selected on the reference field of the filters.
    :param division_value: division selected on the division field of the filters.
    :return: dataset with the filtered data
    """
    filtered_data = unfiltered_data.copy()
    if reference_value == "low + intermediate": reference_value = "rest" #TODO: fix this


     # Apply filters based on dropdown values
    if disease_value:
        filtered_data = filtered_data[filtered_data['score'] == disease_value]
    if quartile_value:
        filtered_data = filtered_data[filtered_data['Quartile'] == quartile_value]
    if reference_value:
        filtered_data = filtered_data[filtered_data['reference'] == reference_value]
    if division_value:
        filtered_data = filtered_data[filtered_data['division'] == division_value]

    return filtered_data

# Callback to change metadata info
@app.callback(
    Output('graph-footer', 'children'),
    [Input('disease-dropdown', 'value'), Input('url', 'pathname')]
)
def update_graph_footer(target_gwas, pathname):
    """
    Callback function to update the graph footer depending on the selected condition.

    :param target_gwas: name of the condition selected.
    :return: formatted text containing the metadata to show on the footer.
    """  

    # Load the result into a Pandas DataFrame
    df = db_handler.get_all_gwas_metadata()
    df = df[df['label'] == target_gwas]

    # Check if we got any result
    if not df.empty:
        # Extract the data from the DataFrame
        row = df.iloc[0] 
        
        # Retrieve date information from metadata
        link_p = row['source']
        link_s = row['sumstats_source']
        population = row['population']
        n = row['n']

        # Create a clickable link
        parts = [f"GWAS information. Sources: "]

        if not pd.isna(link_p):
            parts.append(html.A("Paper", href=link_p, target='_blank'))
        if not pd.isna(link_s):
            if len(parts) > 1:
                parts.append(", ")
            parts.append(html.A("GWAS summary statistics", href=link_s, target='_blank'))

        # Add population and n information
        parts.append(html.Br())  # Add a line break
        parts.append(f"Population: {population}, N: {n}")

        return html.Div(parts)

    else:
        return "No metadata available for the selected GWAS."

# Callback to save clicked point data
@app.callback(
    Output('clicked-data-store', 'data'),
    Input('correlations-graph', 'clickData')
)
def update_clicked_data(clickData):
    if clickData:
        pt = clickData['points'][0]
        # prefer code from customdata (set on the figure); fall back to y
        code = None
        cd = pt.get('customdata')
        if cd is not None:
            # customdata may be scalar or list-like
            try:
                code = cd[0]
            except Exception:
                code = cd
        else:
            code = pt.get('y')
        return {'name': pt.get('x'), 'code': code}
    return no_update

# Callback to update statistics title
@app.callback(
    Output('target-statistics-title', 'children'),
    [Input('clicked-data-store', 'data'),
     State('tabs', 'value')]
)
def update_statistics_title(clicked_data, tab):
    if tab == 'met': return html.H2("Target Statistics", className="section-heading")
    if clicked_data and 'name' in clicked_data:
        return html.H2(f"Target Statistics: {clicked_data['name']}", className="section-heading")
    return html.H2("Target Statistics", className="section-heading")

# Callback to update target information box
@app.callback(
    Output('target-information-box', 'children'),
    [Input('clicked-data-store', 'data'),
     State('tabs', 'value')]
)
def update_basic_statistics(clicked_data, tab):
    # Define a reusable component for the default message
    def default_message():
        return html.Div([
            html.H2("Target Information"),
            html.P("Select a Phecode or ICD code from the graph")
        ], className="info-box")

    # Return the default message if the 'met' tab is selected
    if tab == 'met' or tab == 'quest':
        return default_message()

    # Check if clicked_data contains necessary information
    if not clicked_data or 'name' not in clicked_data:
        return default_message()

    target_desc = clicked_data['code']

    # Execute queries and load results into DataFrames
    prev_target = db_handler.get_disease_prevalence_by_target(target_desc)

    # Ensure data is not empty
    if prev_target.empty:
        return default_message()

    total_individuals = prev_target[prev_target['sex'] == 'both']['prevalence'].values[0]
    total_males = prev_target[prev_target['sex'] == 'male']['prevalence'].values[0]
    total_females = prev_target[prev_target['sex'] == 'female']['prevalence'].values[0]

    # Create a DataFrame for the table
    table_data = pd.DataFrame({
        'Category': ['All', 'Males', 'Females'],
        'Prevalence': [f"{total_individuals:.4f}", f"{total_males:.4f}", f"{total_females:.4f}"]
    })

    # Generate the table
    table = dash_table.DataTable(
        data=table_data.to_dict('records'),
        columns=[{'name': i, 'id': i} for i in table_data.columns],
        style_table={'overflowX': 'auto'},
        style_cell={'textAlign': 'center', 'padding': '10px'},
        style_header={'fontWeight': 'bold'},
        page_size=10
    )
    
    # Add the redirect button
    button = html.Button("View in Cohort", id='redirect-button', n_clicks=0, className='button')

    return html.Div([
        html.H2(f"Statistics for {target_desc}"),
        html.Br(),
        table,
        html.Br(),
        button  # Add the button here
    ], className="info-box")

@app.callback(
    Output('prevalences-graph', 'figure'),
    [
        Input('clicked-data-store', 'data'),
        Input('disease-dropdown', 'value')
    ]
)
def update_prevalences_graph(clicked, disease_value):

    if not clicked or not clicked.get('code'):
        return px.scatter(title='Click a point on the main plot to show prevalence by percentile')

    prs_code = db_handler.get_gwas_code_from_name(disease_value)
    df = db_handler.get_prevalences(prs_code, clicked['code'])
    ttype = db_handler.get_target_type(clicked['code'])

    if df.empty:
        return px.scatter(title='No prevalence data available for selected target')

    # Map PRS column → display label
    df['Sex'] = df['prs_column'].map({
        'PRS_agg': 'All',
        'PRS_male': 'Male',
        'PRS_female': 'Female'
    })

    df = df.dropna(subset=['Sex'])

    if ttype == 'continuous':
        ylabel = 'Mean value'
    else:
        ylabel = 'Prevalence'

    n_groups = df['percentile'].max() +1
    if n_groups == 100:
        xlabel = "Percentile"
    elif n_groups == 10:
        xlabel = "Decile"
    elif n_groups == 4:
        xlabel = "Quartile"
    else:
        xlabel = f"{n_groups} risk score groups"

    fig = px.scatter(
        df,
        x='percentile',
        y='prevalence',
        color='Sex',
        opacity=0.4,
        template='plotly_white'
    )

    # Hide all scatter points by default (show only LOESS lines)
    fig.for_each_trace(
        lambda t: t.update(visible='legendonly')
    )

    # Apply LOESS
    add_loess(fig, df[df['Sex'] == 'All'], 'percentile', 'prevalence', 'All', True)
    add_loess(fig, df[df['Sex'] == 'Female'], 'percentile', 'prevalence', 'Female', 'legendonly')
    add_loess(fig, df[df['Sex'] == 'Male'], 'percentile', 'prevalence', 'Male', 'legendonly')

    fig.update_layout(
        title={'text': f"{(clicked.get('name') if clicked else '')} — {ylabel} by {xlabel}", 'x': 0.5},
        xaxis_title=f'{xlabel} ({disease_value})',
        yaxis_title=ylabel,
        height=PREVALENCE_HEIGHT,
        margin=dict(l=40, r=40, t=60, b=40)
    )

    return fig

def get_target_type(tab):
    if (tab == 'met'):
        return 'Metabolites'
    elif (tab == 'icd'):
        return 'ICD code'
    elif (tab == 'phe'):
        return 'Phecode'
    elif (tab == 'quest'):
        return ('Other Variables')

# Callback to update targets graph in cohort page
@app.callback(
    [Output('phecode-graph', 'figure'),
     Output('icd-code-graph', 'figure')],
    [Input('gender-filter', 'value')]
)
def update_graphs(gender_filter):
    phecode_fig = get_target_plot('Phecodes', gender_filter)
    icd_code_fig = get_target_plot('ICD_codes', gender_filter)
    return phecode_fig, icd_code_fig

# Callback to update the target specific age distribution graph
@app.callback(
    [Output('age-distribution-target-specific-graph', 'figure'),
     Output('bmi-distribution-target-specific-graph', 'figure'),
     Output('hs-distribution-target-specific-graph', 'figure')],
    Input('target-filter', 'value')
)
def update_target_specific_distribution_graph(selected_target):
    age = get_target_specific_dist_graph(selected_target, 'age')
    bmi = get_target_specific_dist_graph(selected_target, 'bmi')
    hs = get_target_specific_hs_graph(selected_target)

    return age, bmi, hs

def get_target_specific_dist_graph(selected_target, output_type):
    if selected_target is None:
        return go.Figure(
            layout=go.Layout(
                xaxis={"visible": False},
                yaxis={"visible": False},
                annotations=[{
                    "text": "Select target to view graph",
                    "xref": "paper",
                    "yref": "paper",
                    "showarrow": False,
                    "font": {
                        "size": 20
                    },
                    "x": 0.5,
                    "y": 0.5,
                    "xanchor": "center",
                    "yanchor": "middle"
                }],
            )
        )
    
    targets_df = db_handler.get_cohort_distribution(selected_target, output_type)

    if output_type == 'age':
        col_name = 'age_group'
        axis_name = 'Age at first diagnosis'
        df_name = 'age'
    elif output_type == 'bmi':
        col_name = 'bmi_group'
        axis_name = 'BMI'
        df_name = 'bmi'

    # Create the bar graph with three bars: Male, Female, and Total
    fig = go.Figure()

    genders = [['MALE', colors[1]], ['FEMALE', colors[2]], ['both', colors[0]]]

    gender_m = {'MALE': 'Male', 'FEMALE': 'Female', 'both': 'Combined'}
    for gender,color in genders:
        gender_data = targets_df[targets_df['sex'] == gender]
        fig.add_trace(go.Bar(
            x=gender_data['category'],
            y=gender_data['count'],
            name=gender_m[gender],
            marker_color=color,
        ))

    title = f'{axis_name} - {selected_target.capitalize()}'
    font_size = get_dynamic_font_size(title)

    # Update layout for better visualization
    fig.update_layout(
        title={
            'text': title,
            'x': 0.5,  # Center the title
            'xanchor': 'center',
            'yanchor': 'top',
            'font': {'size': font_size}
        },
        barmode='group',
        xaxis_title=f'{axis_name}',
        yaxis_title='Nº Individuals',
        xaxis_tickangle=-45,
        margin=dict(t=60, b=20, l=20, r=20),
        paper_bgcolor='rgba(0, 0, 0, 0)',
        plot_bgcolor='rgba(0, 0, 0, 0)',
        xaxis=dict(autorange=True)
    )

    # Hide the legend
    fig.update_layout(showlegend=False)

    return fig

def get_target_specific_hs_graph(selected_target):
    if selected_target is None:
        return go.Figure(
            layout=go.Layout(
                xaxis={"visible": False},
                yaxis={"visible": False},
                annotations=[{
                    "text": "Select target to view graph",
                    "xref": "paper",
                    "yref": "paper",
                    "showarrow": False,
                    "font": {
                        "size": 20
                    },
                    "x": 0.5,
                    "y": 0.5,
                    "xanchor": "center",
                    "yanchor": "middle"
                }],
            )
        )
    
    targets_df = db_handler.get_cohort_distribution(selected_target, 'self_perceived_hs')
    
    hs_order_f = [cat for cat in hs_order if cat in targets_df['category'].values]

    targets_df = (
        targets_df.set_index("category")
        .loc[hs_order_f]
        .reset_index()
    )

    fig = go.Figure()

    # Add bar for males
    if 'MALE' in targets_df['sex'].unique():
        male_df = targets_df[targets_df['sex'] == 'MALE']
        fig.add_trace(go.Bar(
            x=male_df['category'],
            y=male_df['count'],
            name='Males',
            marker_color=colors[1]
        ))

    # Add bar for females
    if 'FEMALE' in targets_df['sex'].unique():
        female_df = targets_df[targets_df['sex'] == 'FEMALE']
        fig.add_trace(go.Bar(
            x=female_df['category'],
            y=female_df['count'],
            name='Females',
            marker_color=colors[2]
        ))
    
    # Add bar for combined data
    if 'both' in targets_df['sex'].unique():
        combined_df = targets_df[targets_df['sex'] == 'both']
        fig.add_trace(go.Bar(
            x=combined_df['category'],
            y=combined_df['count'],
            name='Combined',
            marker_color=colors[0],
            opacity=0.5  # Make combined bars slightly transparent
        ))

    title = f'Self-Perc. HS - {selected_target.capitalize()}'
    font_size = get_dynamic_font_size(title)

    # Update layout for better visualization
    fig.update_layout(
        title={
            'text': title,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top',
            'font': {'size': font_size}
        },
        barmode='group',
        xaxis_title='Self-Perceived HS',
        yaxis_title='Nº Individuals',
        xaxis_tickangle=-45,
        margin=dict(t=60, b=20, l=20, r=20),
        paper_bgcolor='rgba(0, 0, 0, 0)',
        plot_bgcolor='rgba(0, 0, 0, 0)',
        xaxis=dict(autorange=True)
    )

    # Hide the legend
    fig.update_layout(showlegend=False)

    return fig

# Callback to update the target specific age distribution graph
@app.callback(
    Output('gender-distribution-target-specific-graph', 'figure'),
    Input('target-filter', 'value'),
)
def update_target_specific_gender_distribution_graph(selected_target):
    if selected_target is None:
        return go.Figure(
            layout=go.Layout(
                xaxis={"visible": False},
                yaxis={"visible": False},
                annotations=[{
                    "text": "Select target to view graph",
                    "xref": "paper",
                    "yref": "paper",
                    "showarrow": False,
                    "font": {
                        "size": 20
                    },
                    "x": 0.5,
                    "y": 0.5,
                    "xanchor": "center",
                    "yanchor": "middle"
                }]
            )
        )
    
    targets_df = db_handler.get_cohort_distribution(selected_target, 'gender')

    # Filter data based on the selected target
    total_counts = targets_df['count'].sum()
    try:
        female_count = targets_df[targets_df['sex'] == 'FEMALE']['count'].values[0]
    except IndexError as e:
        female_count = 0
    try:
        male_count = targets_df[targets_df['sex'] == 'MALE']['count'].values[0]
    except IndexError as e:
        male_count = 0

    fig = go.Figure(
        data=go.Pie(
            values=[male_count, female_count],
            labels=["Male", "Female"],
            hole=.6,
            direction='clockwise',
            sort=True,
            marker_colors=[colors[1], colors[2]],
            textinfo='label+percent',
            textposition='outside',
            hoverinfo='label+percent'
        )
    )

    title = f'Gender - {selected_target.capitalize()}'
    font_size = get_dynamic_font_size(title)

    fig.update(layout_showlegend=False)
    fig.update_layout(
        title={
            'text': title,
            'x': 0.5,  # Center the title
            'xanchor': 'center',
            'yanchor': 'top',
            'font': {'size': font_size}
        },
        margin=dict(t=55, b=30, l=0, r=0),
    )

    # Hide the legend
    fig.update_layout(showlegend=False)

    return fig

# Callback to change from prs page to cohort page
@app.callback(
    [Output('url', 'pathname', allow_duplicate=True),
     Output('shared-target-data', 'data')],
    [Input('redirect-button', 'n_clicks')],
    [State('clicked-data-store', 'data'),
     State('tabs', 'value')],
    prevent_initial_call=True
)
def redirect_and_set_filter(n_clicks, clicked_data, tab):
    target_type = get_target_type(tab)
    if n_clicks > 0 and clicked_data and 'name' in clicked_data:
        target_desc = clicked_data['name']
        if target_desc:
            return '/cohort', [target_desc, target_type]
    return dash.no_update, dash.no_update

# Callback to update the dropdown value
@app.callback(
    Output('target-filter', 'value'),
    [Input('phecode-graph', 'clickData'),
     Input('icd-code-graph', 'clickData'),
     Input('shared-target-data', 'data')],
)
def update_dropdown_value(phecode_click_data, icd_code_click_data, shared_data):
    ctx = dash.callback_context

    if not ctx.triggered:
        desc = shared_data[0]
        t_type = shared_data[1]
        code = db_handler.get_target_code(desc, t_type)
    else:
        trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]

        if trigger_id == 'phecode-graph' and phecode_click_data:
            code = phecode_click_data['points'][0]['y']
        elif trigger_id == 'icd-code-graph' and icd_code_click_data:
            code = icd_code_click_data['points'][0]['y']
        else:
            return dash.no_update
    return code


# Clientside callback for scrolling
app.clientside_callback(
    """
    function(n_clicks) {
        if(n_clicks > 0) {
            var element = document.getElementById('prs-heading');
            element.scrollIntoView({ behavior: 'smooth' });
        }
        return null;
    }
    """,
    Output('location', 'href'),  # Dummy output; not actually used
    Input('table-button', 'n_clicks')
)

# Callback for the gwas table
@app.callback(
    [Output('url', 'pathname', allow_duplicate=True), Output('selected-gwas', 'data')],
    Input('gwas-table-content', 'active_cell'),
    State('gwas-table-content', 'data'),
    prevent_initial_call=True
)
def on_gwas_row_click(active_cell, table_data):
    if active_cell:
        row = active_cell['row']
        selected_gwas = table_data[row]['name']
        return '/prs', selected_gwas
    return dash.no_update, dash.no_update


# Callback to update the gwas dropdown
@app.callback(
    Output('disease-dropdown', 'value'),
    Input('selected-gwas', 'data')
)
def update_disease_dropdown(selected_gwas):
    if selected_gwas:
        return selected_gwas
    return dash.no_update


#####################################################################################################################
##################################################### App Layout ####################################################
#####################################################################################################################

# Define the desired order for self_perceived_hs
hs_order = ['Very Bad', 'Bad', 'Fair', 'Good', 'Very good', 'DK/NO']

def get_gender_graph():
    
    # Fetch data
    stats = db_handler.get_cohort_distribution('all', 'gender')

    # Calculate percentages
    total_individuals = stats['count'].sum()
    male_count = stats[stats['sex'] == 'MALE']['count'].values[0]
    female_count = stats[stats['sex'] == 'FEMALE']['count'].values[0]
    male_percentage = round((male_count / total_individuals) * 100, 1)
    female_percentage = round((female_count / total_individuals) * 100, 1)

    fig = go.Figure(
        data=go.Pie(
            values=[male_percentage, female_percentage],
            labels=["Male", "Female"],
            hole=.6,
            direction='clockwise',
            sort=True,
            marker_colors=[colors[1], colors[2]],
            textinfo='label+percent',
            textposition='outside',
            hoverinfo='label+percent'
        )
    )

    fig.update_layout(
        margin=dict(t=60, b=0, l=0, r=0),
    )

    fig.update(layout_title_text='Cohort Gender Distribution',
           layout_showlegend=False)

    return fig

def get_distribution_graph(output_type):
    dist = db_handler.get_cohort_distribution('all', output_type)

    fig = go.Figure()

    if output_type == 'self_perceived_hs':
        # Reorder the DataFrame based on the predefined order
        dist = (
            dist.set_index("category")
            .loc[hs_order]
            .reset_index()
        )

    # Add bar for males
    if 'MALE' in dist['sex'].values:
        male_df = dist[dist['sex'] == 'MALE']
        fig.add_trace(go.Bar(
            x=male_df['category'],
            y=male_df['count'],
            name='Males',
            marker_color=colors[1]
        ))

    # Add bar for females
    if 'FEMALE' in dist['sex'].values:
        female_df = dist[dist['sex'] == 'FEMALE']
        fig.add_trace(go.Bar(
            x=female_df['category'],
            y=female_df['count'],
            name='Females',
            marker_color=colors[2]
        ))

    # Add bar for combined data
    if 'both' in dist['sex'].values:
        both_df = dist[dist['sex'] == 'both']
        fig.add_trace(go.Bar(
            x=both_df['category'],
            y=both_df['count'],
            name='Combined',
            marker_color=colors[0],
            opacity=0.5  # Make combined bars slightly transparent
        ))

    if output_type == 'age':
        title = 'Age'
    elif output_type == 'bmi':
        title = 'BMI'
    elif output_type == 'self_perceived_hs':
        title = 'Self-Perceived HS' 

    fig.update_layout(
        title_text=f'Cohort - {title}',
        barmode='group',  # Group bars by category
        xaxis_title=f'{title} Range',
        yaxis_title='Count',
        xaxis_tickangle=-45,
        margin=dict(t=60, b=20, l=20, r=20),
        paper_bgcolor='rgba(0, 0, 0, 0)',  # Transparent background for the paper
        plot_bgcolor='rgba(0, 0, 0, 0)',  # Transparent background for the plot area
    )

    # Hide the legend
    fig.update_layout(showlegend=False)

    return fig

def get_target_plot(target_type, gender_filter='both'):
    df = db_handler.get_disease_prevalence(target_type, gender_filter)

    # Define colors for gender filter
    color_map = {
        'All': [colors[0]],
        'Male': [colors[1]],
        'Female': [colors[2]]
    }
    
    # Determine color based on gender_filter
    colors_local = color_map.get(gender_filter, color_map['All'])

    df = df.sort_values(by='prevalence', ascending=True)

    # Create the bar plot
    fig = px.bar(
        df,
        x='prevalence',
        y='target_code',
        orientation='h',
        title=f"Top 25 {target_type}s by Prevalence ({gender_filter})",
        labels={'prevalence': 'Prevalence', 'target_code':'Code', 'target_description':'Description'},
        color_discrete_sequence=colors_local,  # Apply color sequence
        hover_data={
            'target_code': True,  # Display target ID
            'target_description': True,  # Display target description
            'prevalence': ':.2f'  # Format prevalence with 2 decimal places
        }
    )

    # Set the background to transparent
    fig.update_layout(
        plot_bgcolor='rgba(0, 0, 0, 0)',  # Plot background
        paper_bgcolor='rgba(0, 0, 0, 0)',  # Overall figure background
        height=600  # Set height to 1.5x taller
    )   
    
    return fig

def format_targets_for_dropdown():
    targets = db_handler.get_target_stats()
    targets = targets[targets['target_class'].isin(['Phecodes', 'ICD_codes'])]  # Filter for relevant target classes

    # Remove duplicates based on 'target_id', 'target_description', and 'target_type'
    unique_targets = targets[['target_id', 'target_description', 'target_class']].drop_duplicates()
    unique_targets = unique_targets.sort_values(by='target_id')

    # Create a list of dictionaries for the dropdown options
    dropdown_options = [
        {
            'label': f"{row['target_id']} - {row['target_description'].capitalize()[:25]}{'...' if len(row['target_description']) > 30 else ''} ({row['target_class']})",
            'value': row['target_id']
        }
        for _, row in unique_targets.iterrows()
    ]

    return dropdown_options

def get_dynamic_font_size(title):
    # Set a base font size and reduce it based on the title length
    base_font_size = 24
    max_font_size = 16  # Maximum allowed font size
    length_threshold = 30  # Define a threshold for title length

    # Calculate the font size dynamically
    if len(title) > length_threshold:
        dynamic_size = max(12, base_font_size - (len(title) - length_threshold) * 0.5)
    else:
        dynamic_size = base_font_size

    # Ensure the font size does not exceed the max allowed size
    return min(dynamic_size, max_font_size)


# Navbar
navbar = dbc.Navbar(
    dbc.Container([
        dbc.NavbarBrand(html.Img(src="/assets/logo.png", height="30px"), className="navbar-brand"),
        dbc.NavbarToggler(id="navbar-toggler"),
        dbc.Collapse(
            dbc.Nav([
                dbc.NavItem(dbc.NavLink("PRS Visualisation", href="/", id='prs-link', className="nav-link active")),
                dbc.NavItem(dbc.NavLink("GWAS Sources", href="/gwas", id='gwas-link', className="nav-link")),
                dbc.NavItem(dbc.NavLink("Cohort Overview", href="/cohort", id='cohort-link', className="nav-link")),
                dbc.NavItem(dbc.NavLink("About the Project", href="/about", id='about-link', className="nav-link"))
            ], className="ml-auto", navbar=True),
            id="navbar-collapse",
            navbar=True,
        ),
    ], fluid=True),
    sticky="top",
)



# PRS page layout
prs_page_layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.Div([
                html.H1([
    			html.Img(src="/assets/PolyGenie.png", height="80px", style={"verticalAlign": "middle"})
		]),
                html.P("Unearthing Genetic Links with Polygenic Scores", className="subheading"),
                html.Hr(className="short-sep-line"),
            ], className="prs-info-box"),
            html.Div([
                html.Div([
                    html.H2("Select Filters", className="section-heading"),
                    html.P('Select PRS:', className="dropdown-title"),
                    dcc.Dropdown(
                        id='disease-dropdown',
                        options=disease_options,
                        placeholder='Select phenotype',
                        clearable=True,
                        value=(gwas_names[0] if gwas_names else None),
                        className="dropdown-box"
                    ),
                    html.P('Reference Group (decile/quartile) for Comparison:', className="dropdown-title"),
                    dcc.Dropdown(
                        id='reference-dropdown',
                        options=reference_options,
                        placeholder='Select reference',
                        clearable=False,
                        value="low",
                        className="dropdown-box"
                    ),
                    html.P('PRS Division:', className="dropdown-title"),
                    dcc.Dropdown(
                        id='division-dropdown',
                        options=division_options,
                        placeholder='Select division',
                        clearable=False,
                        value="10",
                        className="dropdown-box"
                    ),
                ]),
                html.Div([
                    html.Button("Check Results Table", id='table-button', className='button'),
                    dcc.Link(
                        html.Button("Check Available GWAS", id='gwas_button', className='button'),
                        href="/gwas"
                    )], className="button-box d-flex justify-content-center gap-3"),
                dcc.Location(id='location', refresh=False),
            ], className="dropdown-box"),
        ], xs=12, sm=12, md=4, lg=4, xl=4, className="left-container"),
            dbc.Col([
                dcc.Tabs(
                    id="tabs",
                    value=default_tab,
                    children=tabs_children
                ),
                dcc.Graph(
                    id='correlations-graph',
                    style={'height': f'{PHEWAS_HEIGHT}px'},
                    config={'displayModeBar': False}
                ),
                html.P(id='graph-footer', className="graph-footer"),
            ],
            xs=12, sm=12, md=8, lg=8, xl=8,
            className="graph-container")
    ], justify='center', className="dashboard-div"),
    html.Br(),
    # Store component to keep track of clicked data point
    dcc.Store(id='clicked-data-store'),
    html.Br(),    
    dbc.Row([
        # First row with basic statistics and one graph
        dbc.Col([
            dcc.Graph(
                id='prevalences-graph',
                className="prevalences-graph",
                style={'height': f'{PREVALENCE_HEIGHT}px'},
                config={'displayModeBar': False}
            )
        ], width=12, md=7, className="full-height custom-col"),  # 70% width for the graph
        dbc.Col([
            html.Div([
                html.H2("Target Information"),
                html.P("Select a Phecode or ICD code from the graph"),
            ], className="info-box")
        ], width=12, md=5, className="target-information-box custom-col", id='target-information-box')  # 30% width for basic stats
    ], className="row-with-margin"),
    html.Br(id="prs-heading"),
    html.H2([
        html.Span("PRS Association Results"),
    ], className="section-heading heading-line"),
    html.Div([
        html.Button("Download Table", id='download-button', className='button ms-auto'),
    ], className="button-box d-flex download-button"),

    # Components to save the table and trigger the download
    dcc.Download(id="download-dataframe-excel"),  
    dcc.Store(id='table-stored-data'), 
        
    dbc.Row([
        html.Div([
            dash_table.DataTable(
                id='prs-table-content',
                columns = [{"name": col, "id": col} for col in col_headers_tokeep],                
                data=[],  
                editable=False,
                sort_action="native",  # Enable sorting
                filter_action='native',
                sort_mode="multi",  # Allow multi-column sorting
                page_action="native",  # Enable pagination
                page_current=0,  # Start from first page
                page_size=20,  # Number of rows per page

                style_cell={
                    'overflow': 'hidden',
                    'textOverflow': 'ellipsis',
                    'whiteSpace': 'normal',
                    'height': 'auto',
                    'maxWidth': '600px'
                },
                style_cell_conditional=[
                    {'if': {'column_id': 'GWAS_code'}, 'width': '80px', 'minWidth': '60px', 'maxWidth': '120px'},
                    {'if': {'column_id': 'Code'}, 'width': '80px', 'minWidth': '60px', 'maxWidth': '120px'},
                    {'if': {'column_id': 'Description'}, 'textAlign': 'left', 'width': '380px', 'minWidth': '200px', 'maxWidth': '600px'},
                    {'if': {'column_id': 'Domain'}, 'textAlign': 'left', 'width': '260px', 'minWidth': '140px', 'maxWidth': '400px'},
                    {'if': {'column_id': 'Beta'}, 'width': '90px', 'minWidth': '70px', 'maxWidth': '120px', 'textAlign': 'right'},
                    {'if': {'column_id': 'OR'}, 'width': '90px', 'minWidth': '70px', 'maxWidth': '120px', 'textAlign': 'right'},
                    {'if': {'column_id': 'CI_5'}, 'width': '90px', 'minWidth': '70px', 'maxWidth': '120px', 'textAlign': 'right'},
                    {'if': {'column_id': 'CI_95'}, 'width': '90px', 'minWidth': '70px', 'maxWidth': '120px', 'textAlign': 'right'},
                    {'if': {'column_id': 'P'}, 'width': '110px', 'minWidth': '80px', 'maxWidth': '140px', 'textAlign': 'right'},
                ],
                
                style_header={
                'backgroundColor': 'rgba(86, 102, 122, .10)',
                'color': 'black',
                'fontWeight': 'bold',
                'textAlign': 'center'
                },

                style_data_conditional=[
                {
                    'if': {'row_index': 'odd'},
                    'backgroundColor': 'rgba(86, 102, 122, .05)',
                },
                {
                    'if': {'state': 'selected'},
                    'backgroundColor': 'rgba(255, 101, 66, .10)',
                    'border': '1px solid #FF6542',
                },
                {
                    'if': {'state': 'active'},
                    'backgroundColor': 'rgba(255, 101, 66, .10)',
                    'border': '1px solid #FF6542',
                }
            ],

            )
        ], id='table-section')
    ]),
    html.Br(),
], fluid=True, style={'margin-top': '50px'})

# Cohort page layout
cohort_page_layout = dbc.Container([
    html.Br(),
    dbc.Row(
        className="cohort-graphs-container",
        children=[
            dbc.Col(
                html.Div([
                    dcc.Graph(
                        figure=get_gender_graph(),
                        className="cohort-graph doghnut-graph",
                    ),
                ]),
                width=12, md=6, lg=3  # Adjust column sizes for responsiveness
            ),
            dbc.Col(
                html.Div([
                    dcc.Graph(
                        figure=get_distribution_graph('age'),
                        id='age-distribution-chart',
                        className="cohort-graph",
                    )
                ]),
                width=12, md=6, lg=3
            ),
            dbc.Col(
                html.Div([
                    dcc.Graph(
                        figure=get_distribution_graph('bmi'),
                        id='bmi-distribution-chart',
                        className="cohort-graph",
                    )
                ]),
                width=12, md=6, lg=3
            ),
            dbc.Col(
                html.Div([
                    dcc.Graph(
                        figure=get_distribution_graph('self_perceived_hs'),
                        id='hs-distribution-chart',
                        className="cohort-graph",
                    )
                ]),
                width=12, md=6, lg=3
            ),
        ]
    ),
    html.Br(),
    html.H2([
        html.Span("Targets Overview"),
    ], className="section-heading heading-line"),
    html.Br(),
    dbc.Row([
        dbc.Col([
            dbc.Row([
                dbc.Col([
                    dcc.Dropdown(
                        id='gender-filter',
                        options=[
                            {'label': 'All', 'value': 'both'},
                            {'label': 'Female', 'value': 'Female'},
                            {'label': 'Male', 'value': 'Male'},
                        ],
                        value='both',
                        placeholder="Filter by Gender",
                        className='dropdown-box cohort-dropdown'
                    ),
                ]),
            ]),
            html.Br(),
            html.Div([
                dcc.Tabs(id="tabs-example", value='tab-1', children=[
                    dcc.Tab(label='Prevalent Phecodes', value='tab-1', children=[
                        dcc.Graph(id='phecode-graph')
                    ]),
                    dcc.Tab(label='Prevalent ICD Codes', value='tab-2', children=[
                        dcc.Graph(id='icd-code-graph')
                    ]),
                ]),
            ], className="phe-prev-bar-graphs"),
        ], width=6, className="custom-col"),
        dbc.Col([
            dbc.Row([
                dbc.Col([
                    dcc.Dropdown(
                        id='target-filter',
                        options=format_targets_for_dropdown(),
                        value='All',
                        placeholder="Filter by Target",
                        className="dropdown-box cohort-dropdown"
                    ),
                ]), 
            ]),
            html.Br(),
            dbc.Row(
                className='right-column-graph-row',
                children=[
                    dbc.Col(
                        html.Div([
                            dcc.Graph(
                                id='age-distribution-target-specific-graph',
                                className="targets-graph",
                            )
                        ]),
                        width=12, md=6  # Full width on small screens, 6 on medium
                    ),
                    dbc.Col(
                        html.Div([
                            dcc.Graph(
                                id='gender-distribution-target-specific-graph',
                                className="targets-graph",
                            )
                        ]),
                        width=12, md=6
                    ),
                ]
            ),
            html.Br(),
            dbc.Row(
                className='right-column-graph-row',
                children=[
                    dbc.Col(
                        html.Div([
                            dcc.Graph(
                                id='bmi-distribution-target-specific-graph',
                                className="targets-graph",
                            )
                        ]),
                        width=12, md=6
                    ),
                    dbc.Col(
                        html.Div([
                            dcc.Graph(
                                id='hs-distribution-target-specific-graph',
                                className="targets-graph",
                            )
                        ]),
                        width=12, md=6
                    ),
                ]
            ),
        ], className="custom-col"),
    ], className='targets-overview-row'),
], fluid=True)

# GWAS page layout

# List of columns to display with custom names
selected_columns = [
    {"name": "GWAS", "id": "name"},
    {"name": "Link Paper", "id": "source"},
    {"name": "Link Summary Statistics", "id": "sumstats_source"},
    {"name": "N", "id": "n"},
]

gwas_page_layout = dbc.Container([
    html.Br(id="prs-heading"),
    html.H2([
        html.Span("PRS Source"),
    ], className="section-heading heading-line"),
    dbc.Row([
        html.Div([
            dash_table.DataTable(
                id='gwas-table-content',
                columns=selected_columns,
                data=db_handler.get_all_gwas_metadata()[[selected_column['id'] for selected_column in selected_columns]].to_dict(orient='records'),  
                editable=False,
                sort_action="native",  # Enable sorting
                filter_action='native',
                sort_mode="multi",  # Allow multi-column sorting
                page_action="native",  # Enable pagination
                page_current=0,  # Start from first page
                page_size=30,  # Number of rows per page
            )
        ], id='table-section')
    ]),
], fluid=True)

# General app layout
app.layout = dbc.Container([
    dcc.Location(id='url', refresh=False),
    navbar,
    # Modal (Popup) for Beta Notice
    html.Div(
        id='terms-modal',
        className='modal',
        children=[
            html.Div(
                className='modal-content',
                children=[
                    html.H4('Usage and Contact Information'),
                    html.P([
                        "For questions, feedback, request a new PRS, or collaboration inquiries, please contact us at: ",
                        html.A("gcatbiobank@igtp.cat", href="mailto:gcatbiobank@igtp.cat")
                    ]),
                    html.P("If the tool contributes to a publication, please cite PolyGenie and the original GWAS used to compute the PRS."),
                    html.Button('Accept', id='accept-button', className='button', n_clicks=0),
                ]
            )
        ]
    ),

    html.Div(id='page-content'),
    dcc.Store(id='shared-target-data', data=[None, None], storage_type='session'),
    dcc.Store(id='selected-gwas', data=None)
], fluid=True)

# Callback to update the page content based on URL
@app.callback(
    Output('page-content', 'children'),
    [Input('url', 'pathname')]
)
def display_page(pathname):
    if pathname == '/cohort':
        return cohort_page_layout
    if pathname == '/about':
        return about_page_layout
    if pathname == '/gwas':
        return gwas_page_layout
    else:
        return prs_page_layout

@app.callback(
    Output('prs-link', 'className'),
    Output('gwas-link', 'className'),
    Output('cohort-link', 'className'),
    Output('about-link', 'className'),
    Input('url', 'pathname')
)
def update_nav_links(pathname):
    if pathname in ['/', '/prs']:
        return 'nav-link active', 'nav-link', 'nav-link', 'nav-link'
    elif pathname == '/gwas':
        return 'nav-link', 'nav-link active', 'nav-link', 'nav-link'
    elif pathname == '/cohort':
        return 'nav-link', 'nav-link', 'nav-link active', 'nav-link'
    elif pathname == '/about':
        return 'nav-link', 'nav-link', 'nav-link', 'nav-link active'

# Callback to hide modal and show content after accepting the terms
@app.callback(
    [Output('terms-modal', 'style'),
     Output('page-content', 'style')],
    [Input('accept-button', 'n_clicks')],
    [State('terms-modal', 'style')]
)
def display_content_after_accept(n_clicks, modal_style):
    if n_clicks > 0:
        # Hide the modal and show the main content
        return {'display': 'none'}, {'display': 'block'}
    # Default state (modal shown, content hidden)
    return {'display': 'block'}, {'display': 'none'}


if __name__ == "__main__":
    app.run_server(debug=True)
