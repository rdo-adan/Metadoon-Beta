import platform
import json
import os
import re
import shutil
import subprocess
from subprocess import CalledProcessError
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import traceback
from PIL import ImageTk
import shutil

#if shutil.which("vsearch") is None:
#    raise EnvironmentError("VSEARCH not found in system PATH. Please install it.")

#if shutil.which("Rscript") is None:
#    raise EnvironmentError("Rscript not found. Please install R.")


class ContainerGeneral:
    metadata = []
    n_threads = os.cpu_count()

class ConfigContainer:
    def __init__(self):
        self.params = {
            "threads": 1,
            "fastq_maxdiffs": 30,
            "fastq_maxee": 1.0,
            "minuniquesize": 2,
            "id": 0.97,
            "sintax_cutoff": 0.8,
            "strand": "both",
            "use_default_chimera_db": True,
            "chimera_db": "",
            "16s_db_option": "rdp",
            "custom_16s_db": ""
        }

    def save_params(self, new_params):
        self.params.update(new_params)
        try:
            save_path = os.path.join(os.getcwd(), 'pipeline_params.json')
            with open(save_path, 'w') as json_file:
                json.dump(self.params, json_file, indent=4)
            messagebox.showinfo("Success", f"Parameters saved to {save_path}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save parameters: {str(e)}")

    def load_params(self):
        try:
            load_path = os.path.join(os.getcwd(), 'pipeline_params.json')
            if os.path.exists(load_path):
                with open(load_path, 'r') as json_file:
                    self.params = json.load(json_file)
                return self.params
            else:
                messagebox.showwarning("Warning", "No parameters file found.")
                return None
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load parameters: {str(e)}")
            return None

config_container = ConfigContainer()

def configure_execution():
    config_window = tk.Toplevel()
    config_window.title("Configure Execution")
    config_window.geometry("400x500")
    config_window.resizable(False, False)

    # Main frame to contain the scrollbar and content
    main_frame = tk.Frame(config_window)
    main_frame.pack(fill="both", expand=True)

    # Canvas to allow scrolling
    canvas = tk.Canvas(main_frame)
    scrollbar = tk.Scrollbar(main_frame, orient="vertical", command=canvas.yview)
    scrollable_frame = tk.Frame(canvas)

    scrollable_frame.bind(
        "<Configure>",
        lambda e: canvas.configure(
            scrollregion=canvas.bbox("all")
        )
    )

    canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
    canvas.configure(yscrollcommand=scrollbar.set)

    def _on_mousewheel(event):
        canvas.yview_scroll(int(-1*(event.delta/120)), "units")
    canvas.bind_all("<MouseWheel>", _on_mousewheel)


    def _on_trackpad(event):
        canvas.yview_scroll(int(-1*event.delta), "units")
    canvas.bind_all("<Shift-MouseWheel>", _on_trackpad)  # Também funciona no Windows
    
    def on_canvas_scroll(event):
        canvas.yview_scroll(int(-1*(event.delta/120)), "units")
    canvas.bind("<MouseWheel>", on_canvas_scroll)


    def on_canvas_focus(event):
        canvas.focus_set()
    canvas.bind("<Enter>", on_canvas_focus)

    # Internal functions
    def toggle_default_chimera_db():
        if default_chimera_db_var.get() == "yes":
            default_chimera_button.config(state=tk.DISABLED)
        else:
            default_chimera_button.config(state=tk.NORMAL)

    def upload_chimera_db():
        file_path_ch = filedialog.askopenfilename(title="Select your Chimera Database")
        if file_path_ch:
            config_container.params["chimera_db"] = file_path_ch

    def select_16s_db():
        option = db_16s_combobox.get()
        if option == "Use Custom 16S Database":
            file_path_16s = filedialog.askopenfilename(title="Select your 16S Database")
            if file_path_16s:
                config_container.params["custom_16s_db"] = file_path_16s
                config_container.params["16s_db_option"] = "custom"
        elif option == "Use RDP 16S Database":
            config_container.params["16s_db_option"] = "rdp"
        else:
            config_container.params["16s_db_option"] = "greengenes2"

    def save_parameters():
        try:
            params = {
                "threads": int(treads.get()),
                "fastq_maxdiffs": int(max_diffs_entry.get()),
                "fastq_maxee": float(max_ee_entry.get()),
                "minuniquesize": int(min_unique_size_entry.get()),
                "id": float(id_entry.get()),
                "sintax_cutoff": float(sintax_cutoff_entry.get()),
                "strand": strand_option.get(),
                "use_default_chimera_db": default_chimera_db_var.get() == "yes",
                "chimera_db": config_container.params.get("chimera_db", ""),
                "metadata": "",
                "tree": "",
                "stat_test": stat_test_combobox.get(),
                "dist_method": dist_method_combobox.get() if 'dist_method_combobox' in locals() else 'bray',
                "color_palette": color_palette_combobox.get(),
                "rarefaction_step": int(rarefaction_step_entry.get()),
                "rarefaction_cex": float(rarefaction_cex_entry.get()),
                "enable_rarefaction": enable_rarefaction_var.get(),
                "rarefaction_depth": int(rarefaction_depth_entry.get()),
                "abundance_top_n": int(abundance_top_n_entry.get()),
                "core_top_n": int(core_top_n_entry.get())
            }
            config_container.save_params(params)
        except ValueError as ve:
            messagebox.showerror("Input Error", f"Please provide valid inputs for all fields. Error: {ve}")

    def load_metadata():
        terminal_output.config(state='normal')
        files = filedialog.askopenfilenames(initialdir="/", title="Select Files")
        for file in files:
            ContainerGeneral.metadata.append(file)
            base_dir = os.getcwd()
            metadata_dir = os.path.join(base_dir, 'Metadata File')
            os.makedirs(metadata_dir, exist_ok=True)
            shutil.copy(file, metadata_dir)
            terminal_output.insert(tk.END, f"Metadata file({file}) loaded...")
        terminal_output.config(state='disabled')
    def load_tree():
        terminal_output.config(state='normal')
        tree_files = filedialog.askopenfilenames(initialdir="/", title="Select Files")
        for tree_file in tree_files:
            ContainerGeneral.metadata.append(tree_file)
            base_dir = os.getcwd()
            tree_dir = os.path.join(base_dir, 'Tree File')
            os.makedirs(tree_dir, exist_ok=True)
            shutil.copy(tree_file, tree_dir)
            terminal_output.insert(tk.END, f"Metadata file({tree_file}) loaded...")
        terminal_output.config(state='disabled')


    # GUI configurations within the scrollable frame
    ttk.Label(scrollable_frame, text="Threads:").pack(pady=5)
    options = list(range(1, ContainerGeneral.n_threads + 1))
    treads = ttk.Combobox(scrollable_frame)
    treads['values'] = options
    treads.set(options[0])
    treads.pack(pady=5)

    ttk.Label(scrollable_frame, text="Max Diffs (Merge):").pack(pady=5)
    max_diffs_entry = ttk.Entry(scrollable_frame)
    max_diffs_entry.insert(0, config_container.params.get("fastq_maxdiffs", 30))
    max_diffs_entry.pack(pady=5)

    ttk.Label(scrollable_frame, text="Max EE (Filter):").pack(pady=5)
    max_ee_entry = ttk.Entry(scrollable_frame)
    max_ee_entry.insert(0, config_container.params.get("fastq_maxee", 1.0))
    max_ee_entry.pack(pady=5)

    ttk.Label(scrollable_frame, text="Min Unique Size (Dereplicate):").pack(pady=5)
    min_unique_size_entry = ttk.Entry(scrollable_frame)
    min_unique_size_entry.insert(0, config_container.params.get("minuniquesize", 2))
    min_unique_size_entry.pack(pady=5)

    ttk.Label(scrollable_frame, text="ID % (Cluster):").pack(pady=5)
    id_entry = ttk.Entry(scrollable_frame)
    id_entry.insert(0, config_container.params.get("id", 0.97))
    id_entry.pack(pady=5)

    ttk.Label(scrollable_frame, text="SINTAX Cutoff:").pack(pady=5)
    sintax_cutoff_entry = ttk.Entry(scrollable_frame)
    sintax_cutoff_entry.insert(0, config_container.params.get("sintax_cutoff", 0.8))
    sintax_cutoff_entry.pack(pady=5)

    ttk.Label(scrollable_frame, text="Strand (Taxonomy):").pack(pady=5)
    strand_option = ttk.Combobox(scrollable_frame, values=["plus", "both"])
    strand_option.set(config_container.params.get("strand", "both"))
    strand_option.pack(pady=5)

    # Chimera database options
    default_chimera_db_var = tk.StringVar(value="yes")
    tk.Checkbutton(scrollable_frame, text="Use Default Chimera Database", variable=default_chimera_db_var, onvalue="yes", offvalue="no", command=toggle_default_chimera_db).pack(pady=5)
    default_chimera_button = tk.Button(scrollable_frame, text="Select Chimera Database", command=upload_chimera_db, state=tk.DISABLED)
    default_chimera_button.pack(pady=5)

    # 16S database options
    ttk.Label(scrollable_frame, text="16S Database Options:").pack(pady=5)
    db_16s_combobox = ttk.Combobox(scrollable_frame, values=["Use RDP 16S Database", "Use Greengenes2 16S Database", "Use Custom 16S Database"])
    db_16s_combobox.set("Use RDP 16S Database")
    db_16s_combobox.pack(pady=5)
    db_16s_combobox.bind("<<ComboboxSelected>>", lambda e: select_16s_db())
    # Upload metadata
    ttk.Label(scrollable_frame, text="Metadata File:").pack(pady=5)
    ttk.Label(scrollable_frame, text="Metadata file here", background="lightgreen", foreground="red").pack(padx=2,pady=2)
    metadata_entry = ttk.Button(scrollable_frame, text="Upload Metadata File", command=load_metadata)
    metadata_entry.pack(pady=30)
    # Upload tree
    ttk.Label(scrollable_frame, text="Tree (Optional):").pack(pady=5)
    ttk.Label(scrollable_frame, text="Tree here", background="lightgreen", foreground="red").pack(padx=2,pady=2)
    tree_entry = ttk.Button(scrollable_frame, text="Upload Tree File", command=load_tree)
    tree_entry.pack(pady=30)
    
    ttk.Label(scrollable_frame, text="Diversity Indices Statistical Test:").pack(pady=5)
    stat_test_combobox = ttk.Combobox(scrollable_frame, values=["t-test", "ANOVA", "Wilcoxon"], state="readonly")
    stat_test_combobox.set(config_container.params.get("stat_test", "ANOVA"))
    stat_test_combobox.pack(pady=5)

    ttk.Label(scrollable_frame, text="Beta Diversity Distance Method:").pack(pady=5)
    dist_method_combobox = ttk.Combobox(scrollable_frame, values=["bray", "jaccard", "euclidean", "manhattan"], state="readonly")
    dist_method_combobox.set(config_container.params.get("dist_method", "bray"))
    dist_method_combobox.pack(pady=5)


    ttk.Label(scrollable_frame, text="Color Palette:").pack(pady=5)
    color_palette_combobox = ttk.Combobox(scrollable_frame, values=["viridis", "plasma", "inferno", "magma", "cividis", "wesanderson"], state="readonly")
    color_palette_combobox.set(config_container.params.get("color_palette", "viridis"))
    color_palette_combobox.pack(pady=5)

    ttk.Label(scrollable_frame, text="Rarefaction Step:").pack(pady=5)
    rarefaction_step_entry = ttk.Entry(scrollable_frame)
    rarefaction_step_entry.insert(0, str(config_container.params.get("rarefaction_step", 100)))
    rarefaction_step_entry.pack(pady=5)

    ttk.Label(scrollable_frame, text="Rarefaction cex:").pack(pady=5)
    rarefaction_cex_entry = ttk.Entry(scrollable_frame)
    rarefaction_cex_entry.insert(0, str(config_container.params.get("rarefaction_cex", 0.6)))
    rarefaction_cex_entry.pack(pady=5)

    #rarefaction
    ttk.Separator(scrollable_frame, orient='horizontal').pack(fill='x', pady=10)
    ttk.Label(scrollable_frame, text="=== Rarefaction Settings ===", font=('Arial', 10, 'bold')).pack(pady=(10,5))
    enable_rarefaction_var = tk.BooleanVar(value=config_container.params.get("enable_rarefaction", False))
    tk.Checkbutton(scrollable_frame, text="Apply rarefaction (subsampling) to normalize data", variable=enable_rarefaction_var).pack(pady=5)
    ttk.Label(scrollable_frame, text="Rarefaction depth (reads per sample):").pack(pady=5)
    rarefaction_depth_entry = ttk.Entry(scrollable_frame, width=15)
    rarefaction_depth_entry.insert(0, str(config_container.params.get("rarefaction_depth", 1000)))
    rarefaction_depth_entry.pack(pady=5)
    ttk.Separator(scrollable_frame, orient='horizontal').pack(fill='x', pady=10)

    ttk.Label(scrollable_frame, text="Rarefaction Depth:").pack(pady=5)
    rarefaction_depth_entry = ttk.Entry(scrollable_frame)
    rarefaction_depth_entry.insert(0, str(config_container.params.get("rarefaction_depth", 1000)))
    rarefaction_depth_entry.pack(pady=5)

    ttk.Label(scrollable_frame, text="Relative Abundance top_n:").pack(pady=5)
    abundance_top_n_entry = ttk.Entry(scrollable_frame)
    abundance_top_n_entry.insert(0, str(config_container.params.get("abundance_top_n", 15)))
    abundance_top_n_entry.pack(pady=5)

    ttk.Label(scrollable_frame, text="Core Microbiome top_n:").pack(pady=5)
    core_top_n_entry = ttk.Entry(scrollable_frame)
    core_top_n_entry.insert(0, str(config_container.params.get("core_top_n", 30)))
    core_top_n_entry.pack(pady=5)
    
    # Button to save the configured parameters
    save_button = ttk.Button(scrollable_frame, text="Save Parameters", command=save_parameters, underline= 100, width=50, takefocus= 125)
    save_button.pack(pady=20)

    # Configure and display the canvas and scrollbar
    canvas.pack(side="left", fill="both", expand=True)
    scrollbar.pack(side="right", fill="y")

# Class to store loaded files
class FileContainer:
    def __init__(self):
        self.files = []
        self.stats = {}

    def add_file(self, file_path):
        if file_path not in self.files:
            self.files.append(file_path)

    def add_stats(self, file_path, stats):
        self.stats[file_path] = stats

container = FileContainer()

class FileContainer:
    def __init__(self):
        self.files = []
        self.stats = {}

    def add_file(self, file_path):
        if file_path not in self.files:
            self.files.append(file_path)

    def add_stats(self, file_path, stats):
        self.stats[file_path] = stats

def generate_stats():
    terminal_output.config(state='normal')
    for file in container.files:
        if file not in container.stats:
            try:
                terminal_output.insert(tk.END, f"Generating statistics for {os.path.basename(file)}...\n")
                command = ["vsearch", "--fastq_eestats2", file, "--output", "-"]
                result = subprocess.run(command, check=True, capture_output=True, text=True)
                container.add_stats(file, result.stdout)
                terminal_output.insert(tk.END, f"Statistics generated for {os.path.basename(file)}\n")
            except FileNotFoundError:
                terminal_output.insert(tk.END, "vsearch not found in the system.\n")
                print("vsearch not found in the system.")
            except CalledProcessError as e:
                terminal_output.insert(tk.END, f"Error calculating statistics for {file}:\n{e.stderr}\n")
                print(f"Error calculating statistics for {file}:\n{e.stderr}\n")
            except Exception as e:
                terminal_output.insert(tk.END, f"Error executing vsearch: {e}\n{traceback.format_exc()}\n")
                print(f"Error executing vsearch: {e}\n{traceback.format_exc()}\n")
    terminal_output.config(state='disabled')

def display_stats():
    stats_text.delete(1.0, tk.END)
    for file in container.files:
        if file in container.stats:
            stats_text.insert(tk.END, f"Statistics for {os.path.basename(file)}:\n")
            stats_text.insert(tk.END, f"{container.stats[file]}\n")
        else:
            stats_text.insert(tk.END, f"Statistics for {file} not found.\n")

def open_file():
    import tkinter as tk
    files = filedialog.askopenfilenames(initialdir="/", title="Select Files",
                                                           filetypes=(("FASTQ Files", "*.fastq"), ("All Files", "*.*")))
    for file in files:
        container.add_file(file)
        file_list.insert(tk.END, file)

    generate_stats()
    display_stats()
    terminal_output.config(state='disabled')

def run_pipeline():
    try:
        try:
            terminal_output.config(state='normal')
            terminal_output.insert(tk.END, "Starting Pipeline...\n")
        except:
            terminal_output.insert(tk.END, " .\n")

        # 1. Merge Pairs - Merging pairs of reads (R1 and R2)
        files_names = []
        merged_name_files = []
        derep_files = []
        clustered_OTUs = []
        # Output directories for the different stages of the pipeline
        base_dir = os.getcwd()  # Use the current directory as base
        merged_dir = os.path.join(base_dir, 'Merged')
        filter_dir = os.path.join(base_dir, 'Filtered')
        derep_dir = os.path.join(base_dir, 'Dereplicated')
        otu_dir = os.path.join(base_dir, 'OTUs')
        tax_dir = os.path.join(base_dir, 'Taxonomy')
        db_dir = os.path.join(base_dir, 'DB')

        # Creating directories as needed
        os.makedirs(merged_dir, exist_ok=True)
        os.makedirs(filter_dir, exist_ok=True)
        os.makedirs(derep_dir, exist_ok=True)
        os.makedirs(otu_dir, exist_ok=True)
        os.makedirs(tax_dir, exist_ok=True)
        os.makedirs(db_dir, exist_ok=True)

        # Load parameters
        params = config_container.load_params()

        # Download Chimera database, if necessary
        if params["use_default_chimera_db"]:
            chimera_db = os.path.join(db_dir, 'rdp_gold.fa')
            if not os.path.exists(chimera_db):
                terminal_output.insert(tk.END, "Downloading default Chimera DB...\n")
                os.system(f"wget http://drive5.com/uchime/rdp_gold.fa --no-check-certificate -O {chimera_db}")
        else:
            chimera_db = params["chimera_db"]

        # Selection of the 16S database
        if params["16s_db_option"] == "rdp":
            sintax_db = os.path.join(db_dir, 'rdp_16s_v16.fa.gz')
            if not os.path.exists(sintax_db):
                terminal_output.insert(tk.END, "Downloading RDP 16S Database...\n")
                os.system(f"wget https://www.drive5.com/sintax/rdp_16s_v16.fa.gz --no-check-certificate -O {sintax_db}")
        elif params["16s_db_option"] == "greengenes2":
            sintax_db = os.path.join(db_dir, '2024.09.seqs.fna.gz')
            if not os.path.exists(sintax_db):
                terminal_output.insert(tk.END, "Downloading Greengenes2 16S Database...\n")
                os.system(f"wget http://ftp.microbio.me/greengenes_release/current/2024.09.seqs.fna.gz --no-check-certificate -O {sintax_db}")
        else:
            sintax_db = params["custom_16s_db"]

        # 1. Merging Reads Forward and Reverse
        for file in container.files:
            if "_R1_" in file:
                file_name = os.path.splitext(os.path.basename(file))[0].replace("_R1_", "_")
                reverse_file = file.replace('_R1_', '_R2_')
                merged_output = os.path.join(merged_dir, f"{file_name}_merged.fastq")
                cmd_merge = f"vsearch --fastq_mergepairs '{file}' --reverse '{reverse_file}' --fastq_maxdiffs 30 --fastq_maxdiffpct 10 --fastqout '{merged_output}' --relabel {file_name}_seq_ --fastq_eeout"
                terminal_output.insert(tk.END, f"\nMerging paired reads for {file_name} and {reverse_file}...\n")
                terminal_output.insert(tk.END, f"Running:{cmd_merge}\n")
                os.system(cmd_merge)
                merged_name_files.append(merged_output)
                files_names.append(file_name)
        for merged in merged_name_files:
            print("Merged Files: \n")
            print(merged)
        # Create directory for combined sample files
        os.makedirs("FullFiles", exist_ok=True)
        subprocess.run(f"cat {merged_dir}/*_merged.fastq > FullFiles/all_single_sample.fastq", shell=True, check=True)
        print("Files combined successfully.")
        # Filtering the reads
        terminal_output.insert(tk.END, "\nFiltering reads...\n")
        cmd_filter = f"vsearch --fastq_filter 'FullFiles/all_single_sample.fastq' --fastq_maxee {params['fastq_maxee']} --sizein --sizeout --fastq_maxns 0 --fastaout {filter_dir}/all_filtered.fasta --fasta_width 0"
        terminal_output.insert(tk.END, f"Running:\n\n {cmd_filter}\n")
        os.system(cmd_filter)
        terminal_output.insert(tk.END, "Reads filtered successfully.\n")

        # Show sequence count at sample level
        print("\nSequences count at sample level:\n")
        os.system(f'cat {filter_dir}/all_filtered.fasta | grep -c "^>"')

        # Further dereplication of all combined samples
        cmd_derep_all = f"vsearch --derep_fulllength '{filter_dir}/all_filtered.fasta' --threads {params['threads']} --sizein --sizeout --output {derep_dir}/all_derep.fasta"
        os.system(cmd_derep_all)

        # Show count of dereplicated sequences
        print("Unique dereplicated sequences: \n")
        os.system(f"grep -c '^>' {derep_dir}/all_derep.fasta")

        # Clustering to create centroids
        cmd_cluster = f"vsearch --cluster_size '{derep_dir}/all_derep.fasta' --threads {params['threads']} --id {params['id']} --strand {params['strand']} --sizein --sizeout --centroids {otu_dir}/samples_centroids.fasta"
        os.system(cmd_cluster)

        # Sort centroids and remove singletons
        cmd_sort = f"vsearch --sortbysize '{otu_dir}/samples_centroids.fasta' --threads {params['threads']} --sizein --sizeout --minsize {params['minuniquesize']} --output {otu_dir}/samples_centroids_sorted.fasta"
        os.system(cmd_sort)

        print("Non-singleton clusters: ")
        os.system(fr'grep -c "^>" {otu_dir}/samples_centroids_sorted.fasta')

        # De novo chimera detection
        cmd_uchime_denovo = f"vsearch --uchime_denovo '{otu_dir}/samples_centroids_sorted.fasta' --sizein --sizeout --fasta_width 0 --qmask none --nonchimeras {otu_dir}/denovo.nonchimeras.fasta"
        os.system(cmd_uchime_denovo)

        print("Unique sequences after de novo chimera detection:")
        os.system(fr'grep -c "^>" {otu_dir}/denovo.nonchimeras.fasta')

        # Reference-based chimera detection
        chimera_db = params["chimera_db"] if not params["use_default_chimera_db"] else os.path.join(db_dir, "rdp_gold.fa")
        cmd_uchime_ref = f"vsearch --uchime_ref '{otu_dir}/denovo.nonchimeras.fasta' --threads {params['threads']} --db '{chimera_db}' --sizein --sizeout --fasta_width 0 --qmask none --dbmask none --nonchimeras {otu_dir}/nonchimeras.fasta"
        os.system(cmd_uchime_ref)

        print("Unique sequences after reference-based chimera detection:")
        os.system(f'grep -c "^>" {otu_dir}/nonchimeras.fasta')

        # Filter and relabel OTUs
        cmd_filter_otus = f"vsearch --fastx_filter '{otu_dir}/nonchimeras.fasta' --threads {params['threads']} --sizein --sizeout --fasta_width 0 --relabel OTU_ --fastaout {otu_dir}/otus.fasta"
        os.system(cmd_filter_otus)

        # Create OTU table
        cmd_otu_table = f"vsearch --usearch_global '{filter_dir}/all_filtered.fasta' --threads {params['threads']} --db '{otu_dir}/otus.fasta' --id {params['id']} --strand {params['strand']} --otutabout {otu_dir}/otutab.txt"
        os.system(cmd_otu_table)

        # Configure SINTAX database
        sintax_db = {
            "rdp": os.path.join(db_dir, "rdp_16s_v16.fa.gz"),
            "greengenes2": os.path.join(db_dir, "2024.09.seqs.fna.gz"),
            "custom": params["custom_16s_db"]
        }[params["16s_db_option"]]

        # Command for taxonomy analysis using SINTAX
        cmd_sintax = f"vsearch --sintax '{otu_dir}/otus.fasta' --db '{sintax_db}' --tabbedout {tax_dir}/taxonomy.txt --sintax_cutoff {params['sintax_cutoff']} --strand {params['strand']}"
        os.system(cmd_sintax)

        # Defines the header to be added
        header = "#OTU ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\t\n"

        # Path to the directory containing the taxonomy.txt file
        file_path = f"{tax_dir}/taxonomy.txt"

        # Check if the file exists
        if os.path.exists(file_path):
            # List to store the processed content
            processed_content = []

            # Read and process the content of the original file
            with open(file_path, 'r') as original_file:
                for line in original_file:
                    # Splits the line into tab-separated columns
                    columns = line.strip().split('\t')

                    # Removes everything after the semicolon in the first column
                    columns[0] = re.sub(r";.*", "", columns[0])

                    # Keeps only the first three columns
                    if len(columns) > 3:
                        columns = columns[:3]

                    # Replaces commas with tabs and removes parentheses and their content in columns after the first
                    columns = [columns[0]] + [re.sub(r"\(.*?\)", "", col.replace(',', '\t').replace('+', '').replace('d:', '').replace('p:', '').replace('c:', '').replace('o:', '').replace('f:', '').replace('g:', '').replace('s:', '')) for col in columns[1:]]

                    # Joins the columns into a modified line and adds it to the processed content list
                    processed_line = '\t'.join(columns) + '\n'
                    processed_content.append(processed_line)

            # Write the new content with the header added
            with open(file_path, 'w') as modified_file:
                modified_file.write(header)
                modified_file.writelines(processed_content)

            print("taxonomy.txt file processed and header added successfully!")

        else:
            print(f"The file {file_path} was not found.")

        try:
            os.system(fr"Rscript Analise.R") # executes the analysis script
            terminal_output.insert(tk.END, "Analysis completed successfully.\n")
        except subprocess.CalledProcessError as e:
            terminal_output.insert(tk.END, f"Error in R analysis execution:\n{e.stderr}\n")  # shows detailed errors
            return  # stops execution if an error occurs
        terminal_output.config(state='normal')
        terminal_output.insert(tk.END, "Pipeline Starting...\n")
    except Exception as e:
        terminal_output.insert(tk.END, f"An unexpected error occurred:\n")
        terminal_output.insert(tk.END, f"{traceback.format_exc()}\n")  # Prints the full traceback
        print(traceback.format_exc())  # Prints the traceback to the console as well
def Generate_report():
    try:
        # Final Report
        result = subprocess.run(
            ["Rscript", "generate_report.R"],
            check=True,
            capture_output=True,
            text=True
        )
        terminal_output.insert(tk.END, "Report Generated Successfully!\n")
        terminal_output.insert(tk.END, result.stdout + "\n")  # opcional: exibe a saída padrão do R
    except subprocess.CalledProcessError as e:
        terminal_output.insert(tk.END, f"Error in Generating the Final Report:\n{e.stderr}\n")
    
    terminal_output.config(state='normal')

def save_analysis_results():
    # Ask user to select the destination folder
    destination = filedialog.askdirectory(title="Select destination folder to save analysis results")

    if not destination:
        terminal_output.insert(tk.END, "Saving canceled by user.\n")
        return

    # Get current working directory (root of the project)
    current_dir = os.getcwd()
    terminal_output.insert(tk.END, f"Saving results from: {current_dir}\n")

    # Get all folders in the current directory
    folders_to_copy = [f for f in os.listdir(current_dir) if os.path.isdir(f)]

    # Copy each folder to the selected destination
    for folder in folders_to_copy:
        src = os.path.join(current_dir, folder)
        dst = os.path.join(destination, folder)
        try:
            if os.path.exists(dst):
                shutil.rmtree(dst)  # Remove if already exists
            shutil.copytree(src, dst)
            terminal_output.insert(tk.END, f"Copied folder: {folder}\n")
        except Exception as e:
            terminal_output.insert(tk.END, f"Error copying {folder}: {e}\n")

    # Copy specific important files (if they exist)
    for file_name in ["Rplots.pdf", "pipeline_params.json", "Metadoon_Report.html"]:
        file_path = os.path.join(current_dir, file_name)
        if os.path.exists(file_path):
            try:
                shutil.copy(file_path, destination)
                terminal_output.insert(tk.END, f"Copied file: {file_name}\n")
            except Exception as e:
                terminal_output.insert(tk.END, f"Error copying {file_name}: {e}\n")

    # ✅ Cleanup: remove original folders and files only if saving was successful
    for folder in folders_to_copy:
        src = os.path.join(current_dir, folder)
        try:
            shutil.rmtree(src)
            terminal_output.insert(tk.END, f"Removed original folder: {folder}\n")
        except Exception as e:
            terminal_output.insert(tk.END, f"Error removing folder {folder}: {e}\n")

    for file_name in ["Rplots.pdf", "pipeline_params.json", "Metadoon_Report.html"]:
        file_path = os.path.join(current_dir, file_name)
        if os.path.exists(file_path):
            try:
                os.remove(file_path)
                terminal_output.insert(tk.END, f"Removed original file: {file_name}\n")
            except Exception as e:
                terminal_output.insert(tk.END, f"Error removing file {file_name}: {e}\n")

    terminal_output.insert(tk.END, "✅ Results copied and original files cleaned up successfully.\n")


# MAIN WINDOW
root = tk.Tk()
root.title("Metadoon Version 1.0")

# 🚀 Detecta sistema operacional e define o ícone adequado
system_os = platform.system()

if system_os == "Windows":
    root.iconbitmap("Metadoon.ico")
elif system_os == "Darwin":  # macOS
    img = ImageTk.PhotoImage(file="Metadoon.icns")
    root.iconphoto(True, img)
else:  # Linux e outros
    img = ImageTk.PhotoImage(file="Metadoon.png")
    root.iconphoto(True, img)

# ✅ Continua seu código normalmente
container = FileContainer()

# Menus
menu = tk.Menu(root)
root.config(menu=menu)

file_menu = tk.Menu(menu, tearoff=0)
menu.add_cascade(label="Files", menu=file_menu)
file_menu.add_command(label="Upload Files", command=open_file)
file_menu.add_command(label="Exit", command=root.quit)

tools_menu = tk.Menu(menu, tearoff=0)
menu.add_cascade(label="Tools", menu=tools_menu)
tools_menu.add_command(label="Configure Execution", command=configure_execution)
tools_menu.add_command(label="Run Pipeline", command=run_pipeline)
tools_menu.add_command(label="Generate Final Report", command=Generate_report)
tools_menu.add_command(label="Save and Clean Results", command=save_analysis_results)

# Dividing the window into two columns (Files and Statistics)
file_label = tk.Label(root, text="FILES", bg="lightgreen", foreground="black")
file_label.grid(row=0, column=0, padx=10, pady=5, sticky="nsew")

file_list = tk.Listbox(root)
file_list.grid(row=1, column=0, padx=10, pady=5, sticky="nsew")

stats_label = tk.Label(root, text="FILE STATS:", bg="lightgreen", foreground="black")
stats_label.grid(row=0, column=1, padx=10, pady=5, sticky="nsew")

stats_text = tk.Text(root)
stats_text.grid(row=1, column=1, padx=10, pady=5, sticky="nsew")

terminal_label = tk.Label(root, text="STATUS", bg="lightgreen", foreground="red")
terminal_label.grid(row=2, column=0, columnspan=2, padx=10, pady=5, sticky="nsew")

terminal_output = tk.Text(root, height=5)
terminal_output.grid(row=3, column=0, columnspan=2, padx=10, pady=5, sticky="nsew")

progress_bar = ttk.Progressbar(root)
progress_bar.grid(row=4, column=0, columnspan=2, padx=10, pady=5, sticky="nsew")

configure_button = tk.Button(root, text="Configure Execution", command=configure_execution, foreground="black")
configure_button.grid(row=5, column=0, padx=10, pady=5, sticky="nsew")

run_button = tk.Button(root, text="RUN PIPELINE", command=run_pipeline, foreground="black")
run_button.grid(row=5, column=1, padx=10, pady=5, sticky="nsew")

root.grid_columnconfigure(0, weight=2)
root.grid_columnconfigure(1, weight=2)
root.grid_rowconfigure(1, weight=1)

root.mainloop()
