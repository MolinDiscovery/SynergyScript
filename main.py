import pandas as pd
import numpy as np
import ast
import os
import re
import glob

from bokeh.models import FactorRange, Span, CustomJS, Switch, Whisker, Button, BuiltinIcon
from bokeh.plotting import figure, show, output_notebook
from bokeh.models.sources import ColumnDataSource
from bokeh.layouts import column
from bokeh.io import push_notebook
import bokeh.palettes as bpcolors
from bokeh.transform import factor_cmap

from bokeh.models import CustomJS, Dropdown, TabPanel, Tabs
from bokeh.plotting import figure
from bokeh.layouts import column, row, Spacer
from bokeh.io import show, curdoc, output_notebook, reset_output, output_file


script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)
directory_path = "Data"
directory_path = os.path.join(script_directory, directory_path)



def read_and_process_files(base_directory):
    """
    Reads files from multiple batch folders and organizes the data into a nested structure.

    This function processes .xlsx and .txt files from batch folders within the specified base directory. 
    It organizes the data into a nested structure with each batch folder as a key.

    Parameters:
    base_directory (str): The path to the base directory containing batch folders.

    Returns:
    dict: A nested dictionary with each batch folder as a key and its data structure as the value.
    """

    def process_xlsx(file):
        """
        Reads an Excel file and processes the data into a DataFrame.
        """
        df = pd.read_excel(file, index_col=0)
        df = df.iloc[:8, :12]  # Keep only the first 8x12 rows x columns
        df = df.applymap(lambda x: x / 1000 if x >= 10 else x)
        return df

    def process_txt(file):
        """
        Reads a text file and extracts concentration data.
        """
        with open(file, 'r') as f:
            lines = f.readlines()
        conc_lists = [ast.literal_eval(line.strip()) for line in lines if line.strip()]
        return conc_lists

    all_batches_data = {}

    # Iterate over batch folders
    batch_folders = glob.glob(os.path.join(base_directory, '*_batch'))
    for batch_folder in batch_folders:
        synergy_folder = os.path.join(batch_folder, 'synergy')

        # Initialize the data structure for this batch
        batch_data = {'data': [], 'conc': []}
        data_frames = {}

        # Process each .xlsx file in the synergy folder
        for xlsx_file in glob.glob(os.path.join(synergy_folder, '*.xlsx')):
            match = re.match(r'(\d+) - (\d+) - Raw Data_(\d+)\.xlsx', os.path.basename(xlsx_file))
            if match:
                exp_num, rep_num, _ = match.groups()
                if (exp_num, rep_num) not in data_frames:
                    data_frames[(exp_num, rep_num)] = []
                data_frames[(exp_num, rep_num)].append(process_xlsx(xlsx_file))

        # Process each .txt file in the synergy folder
        for txt_file in glob.glob(os.path.join(synergy_folder, '*.txt')):
            match = re.match(r'(\d+) - (\d+) - conc\.txt', os.path.basename(txt_file))
            if match:
                exp_num, rep_num = match.groups()
                conc_data = process_txt(txt_file)
                batch_data['conc'].append((int(exp_num), int(rep_num), *conc_data))
                tuples = [(int(rep_num), int(exp_num), df) for df in data_frames.get((exp_num, rep_num), [])]
                batch_data['data'].append(tuples)

        # Assign the batch data to the all_batches_data dictionary
        batch_name = os.path.basename(batch_folder)
        all_batches_data[batch_name] = batch_data

    return all_batches_data


def read_data(data_files: list = glob.glob('*.xlsx')):

    """
    Reads all .xlsx files from the provided list of files and returns a list of dataframes.
    Only the first 8x12 rows x columns of each dataframe are kept.

    Parameters:
    data_files (list): A list of .xlsx file paths. Defaults to all .xlsx files in the current directory.
    
    Returns:
    list: A list of pandas DataFrames, each DataFrame corresponding to an Excel file.
    """

    df_list = []

    for file in data_files:
        df = pd.read_excel(file, index_col=0)
        df = df.iloc[:8, :12] # Only keep the first 8x12 rows x columns
        df = df.applymap(lambda x: x / 1000 if x >= 10 else x)
        df_list.append(df)

    return df_list

def format_to_synergy(catagorized_data: pd.DataFrame, Drug1_name = 'Drug 1', Drug2_name = 'Drug 2'):
    
    df = catagorized_data

    def extract_concentrations(treatment_str):
        """
        Extracts drug concentrations from the treatment string.
        Assumes the format 'd1d2 (x / y)' where x and y are concentrations.
        """
        if 'd1d2' in treatment_str:
            parts = treatment_str.split('(')[-1].split(')')[0].split('/')
            drug1_conc = float(parts[0].strip())
            drug2_conc = float(parts[1].strip())
        elif 'd1' in treatment_str:
            parts = treatment_str.split('(')[-1].split(')')[0]
            drug1_conc = float(parts.strip())
            drug2_conc = 0
        elif 'd2' in treatment_str:
            parts = treatment_str.split('(')[-1].split(')')[0]
            drug1_conc = 0
            drug2_conc = float(parts.strip())
        else:
            # Handle other cases or set default values
            drug1_conc = 0
            drug2_conc = 0

        return drug1_conc, drug2_conc
    

    df['Conc1'], df['Conc2'] = zip(*df['treatment'].apply(extract_concentrations))
    df['Response'] = df['control_pct']

    # remove rows where groups column is media(backup), DMSO, BG.
    df = df[~df['groups'].isin(['media(backup)', 'DMSO', 'BG'])]
    df = df[['Conc1', 'Conc2', 'Response']]
    
    df['Drug1'] = Drug1_name
    df['Drug2'] = Drug2_name
    df['PairIndex'] = 1
    df['ConcUnit'] = 'uM'

    df = df[['PairIndex', 'Drug1', 'Drug2', 'Conc1', 'Conc2', 'Response', 'ConcUnit']]
    return df

class ReadPlates:

    def __init__(
            self,
            data_frames: list,
            drug_concentrations: list,
            ) -> None:
        
        self.data_frames = data_frames
        self.num_experiments = len(data_frames)
        self.drug_concentrations = drug_concentrations
        self.plates = None
        self.plates_calculated = None
        self.plates_replicates = None


    def read_plates(self) -> list:

        plates = []
        count = 0

        for df in self.data_frames:
            first_three_cols = df.iloc[:, :3]
            second_three_cols = df.iloc[:, 3:6]
            third_three_cols = df.iloc[:, 6:9]
            second_three_cols.columns = first_three_cols.columns
            third_three_cols.columns = first_three_cols.columns

            drug2_1 = list(df.loc[['A', 'B', 'C'],10])
            drug2_2 = list(df.loc[['D', 'E', 'F'],10])
            drug2_3 = list(df.loc[['A', 'B', 'C'],11])
            drug2_4 = list(df.loc[['D', 'E', 'F'],11])
            drug_df = pd.DataFrame(columns=[1, 2, 3])
            drug_df.loc[1] = drug2_1
            drug_df.loc[2] = drug2_2
            drug_df.loc[3] = drug2_3
            drug_df.loc[4] = drug2_4

            extracted_data = pd.concat([first_three_cols, second_three_cols, third_three_cols, drug_df], axis=0)

            index_list = []
            # items = list(self.drug_concentrations.items())
            # key, value = items[count]

            value = self.drug_concentrations[count]

            for drug2_conc in value[1]:
                for drug1_conc in value[0]:
                    index_list.append(f'd1d2 ({drug1_conc} / {drug2_conc})')

            for drug1_conc in value[0]:
                index_list.append(f'd1({drug1_conc})')

            index_list.append('media (backup)')
            index_list.append('DMSO')
            index_list.append('media')
            index_list.append('BG')

            for drug2_conc in value[1]:
                index_list.append(f'd2({drug2_conc})')
            
            extracted_data.index = index_list

            plates.append(extracted_data)
            count += 1

        self.plates = plates

        return self.plates
    

    def calculate_statistics(self) -> list:

        plates_calculated = []

        if self.plates is None:
            self.read_plates()

        for plate in self.plates:
            plate['avg'] = plate.mean(axis=1)
            plate['std'] = plate.std(axis=1, ddof=1)
            plate['control_pct'] = plate['avg'] / plate.loc['media', 'avg'] * 100
            plate['std_pct'] = plate['std'] / plate.loc['media', 'avg'] * 100

            plates_calculated.append(plate)

        self.plates_calculated = plates_calculated

        return self.plates_calculated
    

    def categorize_data(self, replicates=False) -> list:

        # Calculate the statistics if they haven't been calculated yet
        if self.plates_calculated is None:
            self.plates_calculated = self.calculate_statistics()

        # Calculate the replicates if they haven't been calculated yet
        if self.plates_replicates is None and replicates:
            self.plates_replicates = self.calculate_replicates()
        
        # Create a list to hold the categorized plates
        plates_categorized = []

        
        if replicates:
            # If replicates are calculated, use those plates
            plates = self.plates_replicates
        else:
            # Otherwise, use the calculated plates
            plates = self.plates_calculated

        count = 0
        # Categorize the plates
        for plate in plates:
            plate.loc[plate.index[0:4], 'groups'] = f"Drug 2: {self.drug_concentrations[count][1][0]}"
            plate.loc[plate.index[4:8], 'groups'] = f"Drug 2: {self.drug_concentrations[count][1][1]}"
            plate.loc[plate.index[8:12], 'groups'] = f"Drug 2: {self.drug_concentrations[count][1][2]}"
            plate.loc[plate.index[12:16], 'groups'] = f"Drug 2: {self.drug_concentrations[count][1][3]}"
            plate.loc[plate.index[16:20], 'groups'] = "Drug 1"
            plate.loc[plate.index[24:25], 'groups'] = f"Drug 2: {self.drug_concentrations[count][1][0]}"
            plate.loc[plate.index[25:26], 'groups'] = f"Drug 2: {self.drug_concentrations[count][1][1]}"
            plate.loc[plate.index[26:27], 'groups'] = f"Drug 2: {self.drug_concentrations[count][1][2]}"
            plate.loc[plate.index[27:28], 'groups'] = f"Drug 2: {self.drug_concentrations[count][1][3]}"
            plate.loc[plate.index[20:21], 'groups'] = "media(backup)"
            plate.loc[plate.index[21:22], 'groups'] = "DMSO"
            plate.loc[plate.index[22:23], 'groups'] = "media"
            plate.loc[plate.index[23:24], 'groups'] = "BG"
            plate.reset_index(inplace=True)
            plate.rename(columns={'index':'treatment'}, inplace=True)
            plate.insert(0, 'groups', plate.pop('groups'))
            plate.columns = [str(col) for col in plate.columns]
            
            count += 1

            plates_categorized.append(plate)
        
        return plates_categorized
    

    def calculate_replicates(self) -> list:
            
            if self.plates_calculated is None:
                self.plates_calculated = self.calculate_statistics()
    
            plates_replicates = []
            replicate_plates = pd.DataFrame()
            # extract the contorl_pct column from each plate
            # put these into a new dataframe
            count = 0
            for plate in self.plates_calculated:
                replicate_plates[f'control_pct_{str(count)}'] = list(plate['control_pct'])
                #print(plate['control_pct'])
                count += 1
            
            replicate_plates.index = self.plates_calculated[0].index

            # calculate the mean of each row
            replicate_plates['control_pct'] = replicate_plates.mean(axis=1)
            # calculate the standard error of each row
            replicate_plates['std'] = replicate_plates.std(axis=1, ddof=1)
            replicate_plates['sem'] = replicate_plates['std'] / np.sqrt(replicate_plates.count(axis=1))
    
            plates_replicates.append(replicate_plates)
    
            self.plates_replicates = plates_replicates
    
            return self.plates_replicates

class PlotSynergy:

    def __init__(self, catagorized_data: list) -> None:

        self.catagorized_data = catagorized_data

    def plot_synergy(self, replicates=False, show_plot=False) -> None:
        
        if replicates:
            control_pct = 'control_pct'
            std_pct = 'sem'
        else:
            control_pct = 'control_pct'
            std_pct = 'std_pct'

        all_layouts = []

        for plate in self.catagorized_data:
            media = plate['groups'][22]
            media = plate['groups'][22]
            DMSO = plate['groups'][21]
            drug_1 = plate['groups'][16]
            drug_1_sub = list(plate['treatment'][16:20])

            drug_2_Q1 = plate['groups'][0]
            drug_2_Q1_alone = plate['treatment'][24]
            drug_2_Q1_sub = list(plate['treatment'][0:4])

            drug_2_Q2 = plate['groups'][4]
            drug_2_Q2_alone = plate['treatment'][25]
            drug_2_Q2_sub = list(plate['treatment'][4:8])

            drug_2_Q3 = plate['groups'][8]
            drug_2_Q3_alone = plate['treatment'][26]
            drug_2_Q3_sub = list(plate['treatment'][8:12])

            drug_2_Q4 = plate['groups'][12]
            drug_2_Q4_alone = plate['treatment'][27]
            drug_2_Q4_sub = list(plate['treatment'][12:16])

            data = {
                'groups': [
                    (media, media),
                    (DMSO, DMSO),
                    (drug_1, drug_1_sub[0]),
                    (drug_1, drug_1_sub[1]),
                    (drug_1, drug_1_sub[2]),
                    (drug_1, drug_1_sub[3]),
                    (drug_2_Q1, drug_2_Q1_alone),
                    (drug_2_Q1, drug_2_Q1_sub[0]),
                    (drug_2_Q1, drug_2_Q1_sub[1]),
                    (drug_2_Q1, drug_2_Q1_sub[2]),
                    (drug_2_Q1, drug_2_Q1_sub[3]),
                    (drug_2_Q2, drug_2_Q2_alone),
                    (drug_2_Q2, drug_2_Q2_sub[0]),
                    (drug_2_Q2, drug_2_Q2_sub[1]),
                    (drug_2_Q2, drug_2_Q2_sub[2]),
                    (drug_2_Q2, drug_2_Q2_sub[3]),
                    (drug_2_Q3, drug_2_Q3_alone),
                    (drug_2_Q3, drug_2_Q3_sub[0]),
                    (drug_2_Q3, drug_2_Q3_sub[1]),
                    (drug_2_Q3, drug_2_Q3_sub[2]),
                    (drug_2_Q3, drug_2_Q3_sub[3]),
                    (drug_2_Q4, drug_2_Q4_alone),
                    (drug_2_Q4, drug_2_Q4_sub[0]),
                    (drug_2_Q4, drug_2_Q4_sub[1]),
                    (drug_2_Q4, drug_2_Q4_sub[2]),
                    (drug_2_Q4, drug_2_Q4_sub[3]),
                ],
                'control_pct': (
                    plate[control_pct][22],
                    plate[control_pct][21],
                    plate[control_pct][16],
                    plate[control_pct][17],
                    plate[control_pct][18],
                    plate[control_pct][19],
                    plate[control_pct][24],
                    plate[control_pct][0],
                    plate[control_pct][1],
                    plate[control_pct][2],
                    plate[control_pct][3],
                    plate[control_pct][25],
                    plate[control_pct][4],
                    plate[control_pct][5],
                    plate[control_pct][6],
                    plate[control_pct][7],
                    plate[control_pct][26],
                    plate[control_pct][8],
                    plate[control_pct][9],
                    plate[control_pct][10],
                    plate[control_pct][11],
                    plate[control_pct][27],
                    plate[control_pct][12],
                    plate[control_pct][13],
                    plate[control_pct][14],
                    plate[control_pct][15],
                ),
                'std_pct': (
                    plate[std_pct][22],
                    plate[std_pct][21],
                    plate[std_pct][16],
                    plate[std_pct][17],
                    plate[std_pct][18],
                    plate[std_pct][19],
                    plate[std_pct][24],
                    plate[std_pct][0],
                    plate[std_pct][1],
                    plate[std_pct][2],
                    plate[std_pct][3],
                    plate[std_pct][25],
                    plate[std_pct][4],
                    plate[std_pct][5],
                    plate[std_pct][6],
                    plate[std_pct][7],
                    plate[std_pct][26],
                    plate[std_pct][8],
                    plate[std_pct][9],
                    plate[std_pct][10],
                    plate[std_pct][11],
                    plate[std_pct][27],
                    plate[std_pct][12],
                    plate[std_pct][13],
                    plate[std_pct][14],
                    plate[std_pct][15],
                ),                
            }
            # Extract the second elements of the tuples
            treatment_groups = [group[1] for group in data['groups']]

            # Add the second elements to the data source
            data['treatment_groups'] = treatment_groups            

            # Calculate error bar ends
            data['upper'] = [control + std for control, std in zip(data['control_pct'], data['std_pct'])]
            data['lower'] = [control - std for control, std in zip(data['control_pct'], data['std_pct'])]

            source = ColumnDataSource(data=data)

            unique_groups = list(dict.fromkeys(x[0] for x in data['groups']))
            color = bpcolors.colorblind['Colorblind'][len(unique_groups)]
            index_cmap = factor_cmap('groups', palette=color, factors=unique_groups, end=1)

            tooltip_error_metric = 'Standard error' if replicates else 'Standard deviation'
            p = figure(
                x_range=FactorRange(*data['groups']),
                height=350, width=1000,
                toolbar_location='right',
                tooltips=[('Treatment', '@treatment_groups'),
                          ('Percentage of control', '@control_pct'),
                          (tooltip_error_metric, '@std_pct')],
                title='Synergy Plot')

            p.vbar(x='groups', top='control_pct', width=1, source=source, line_color="white", fill_color=index_cmap,)

            # Add error bars
            p.segment(x0='groups', y0='lower', x1='groups', y1='upper', line_color='black', source=source)

            # Add caps to error bars
            whisker = Whisker(source=source, base='groups', upper='upper', lower='lower', level='overlay', line_width=0.1)
            p.add_layout(whisker)

            span = Span(location=100, dimension='width', line_color='red', line_width=1)
            p.renderers.extend([span])

            ''' Disabled for now
            # Create a Switch widget
            switch = Switch(active=True)

            # Define a CustomJS callback
            callback = CustomJS(args=dict(span=span), code="""
                span.visible = cb_obj.active;
            """)

            # Attach the callback to the Switch widget
            switch.js_on_change('active', callback)
            '''

            '''Creates a downloadable csv file for SynergyFinder'''
            synergy_data = format_to_synergy(plate)
            synergy_data_csv = synergy_data.to_csv(index=False)
            csv_data = synergy_data_csv.encode().decode('utf-8')
            
            icon = BuiltinIcon("settings", size="1.2em", color="white")
            download_button = Button(label='Download Synergy Data', button_type='success', icon=icon)

            # CustomJS to download CSV file
            download_callback = CustomJS(args=dict(csv_data=csv_data), code="""
                var blob = new Blob([csv_data], { type: 'text/csv;charset=utf-8;' });
                var link = document.createElement("a");
                link.href = URL.createObjectURL(blob);
                link.download = "synergy_data.csv";
                link.click();
                URL.revokeObjectURL(link.href);
            """)

            download_button.js_on_click(download_callback)


            p.y_range.start = 0
            p.x_range.range_padding = 0.05
            p.xgrid.grid_line_color = None
            #p.xaxis.axis_label = "Drug are Wonderful"
            p.xaxis.major_label_orientation = 1.3
            p.xaxis.major_label_text_font_size = '0pt' # turn off x-axis tick labels
            p.outline_line_color = None
            p.yaxis.axis_label = "%Viability relative to CTRL"

            # Add the widget to the layout
            #layout = column(p, switch)
            layout = column(p, download_button)
            all_layouts.append(layout)

        if show_plot:
            output_notebook()
            for layout in all_layouts:
                show(layout)  # Display each plot layout

        # Return the list of layouts if not displaying
        if not show_plot:
            return all_layouts

class DashBoard:

    pass

batches = read_and_process_files(directory_path)
#print(batches)

# Initialize variables to hold batch names and compounds in each batch
batch_names = []
compound_in_batch = {}
unique_compounds = []
plate_data = []
# Iterate over batches and extract the compound names
all_plots = {}
for batch_key, batch_values in batches.items():
    batch_names.append((batch_key, batch_key))

    data = batch_values.get('data')
    conc = batch_values.get('conc')

    # Retrieve dataframes for each replicate
    compound_dataframes = []
    plots = []
    count = 0
    for compound in data:
        # Check if there are replicates
        replicate = True if len(compound) > 1 else False

        for tup in compound:
            compound_try = tup[0]
            compound_name = tup[1]
            compound_dataframes.append(tup[2])

        conc_d1 = tuple(conc[count][2])
        conc_d2 = tuple(conc[count][3])
        conccentrations = [conc_d1, conc_d2]
        conccentrations = [conccentrations for _ in range(len(compound_dataframes))] # Repeat the concentrations for each replicate
        count = count + 1
        calculate_stats = ReadPlates(compound_dataframes, conccentrations).categorize_data()
        plot = PlotSynergy(calculate_stats).plot_synergy()
        if replicate:
            replicate_plates = ReadPlates(compound_dataframes, conccentrations).categorize_data(replicates=True)
            replicate_plot = PlotSynergy(replicate_plates).plot_synergy(replicates=True)
            plot.insert(0, replicate_plot[0])
        plots.append((compound_try, compound_name, plot))
        compound_dataframes = []
        all_plots[batch_key] = plots
        replicate = False

    plots = []
    

### Creating dashboard ###    
menu = [(batch, batch) for batch in all_plots.keys()]
dropdown = Dropdown(label="Select Batch", button_type="default", menu=menu)

def create_dynamic_plot_layout(plots, total_width):
    if len(plots) == 0:
        return None

    # Set the large plot
    large_plot = plots[0]
    large_plot.width = total_width
    large_plot.height = 400  # This just adds padding between large and small plots

    layout = [large_plot]  # Start with the large plot
    
    # Process the smaller plots
    small_plots = plots[1:]

    # Create rows of smaller plots
    for i in range(0, len(small_plots), 2):
        row_plots = small_plots[i:i+2]
        layout.append(row(*row_plots))

    return column(*layout)

# Use this function in the construction of tabs_dict with a specified total width
total_layout_width = 4000  # Example total width; adjust as needed
tabs_dict = {
    batch: [
        TabPanel(child=create_dynamic_plot_layout(value[2], total_layout_width), title=f"{value[1]}({value[0]})")
        for value in values
    ]
    for batch, values in all_plots.items()
}

# Use this function in the construction of tabs_dict
initial_batch = list(all_plots.keys())[0]
tabs = Tabs(tabs=tabs_dict[initial_batch])

# CustomJS callback for dropdown
callback = CustomJS(args=dict(tabs=tabs, tabs_dict=tabs_dict), code="""
    const batch = cb_obj.item;
    tabs.tabs = tabs_dict[batch];
""")
dropdown.js_on_event("menu_item_click", callback)

# Add to the current document
layout = column(dropdown, tabs)
curdoc().add_root(layout)

output_file('index.html', title='SynergyScript | MolinDiscovery')

show(layout)