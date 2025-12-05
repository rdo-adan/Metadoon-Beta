import platform
import glob
import json
import os
import re
import shutil
import subprocess
import threading
import time
from subprocess import CalledProcessError
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import traceback
from PIL import ImageTk
from datetime import datetime

# --- DEPENDENCY CHECK ---
def check_dependencies():
    missing = []
    if shutil.which("vsearch") is None:
        missing.append("VSEARCH")
    if shutil.which("Rscript") is None:
        missing.append("R (Rscript)")
    
    if missing:
        messagebox.showwarning("Dependency Error", f"Missing dependencies: {', '.join(missing)}.\nPlease install them and add to PATH.")

class ContainerGeneral:
    metadata = []
    n_threads = os.cpu_count()

# --- LOADING WINDOW ---
class LoadingWindow:
    def __init__(self, parent, title="Processing"):
        self.top = tk.Toplevel(parent)
        self.top.title(title)
        self.top.geometry("350x120")
        self.top.resizable(False, False)
        
        x = parent.winfo_x() + (parent.winfo_width() // 2) - 175
        y = parent.winfo_y() + (parent.winfo_height() // 2) - 60
        self.top.geometry(f"+{x}+{y}")

        self.top.transient(parent)
        self.top.grab_set() 

        tk.Label(self.top, text="Please wait...", font=("Arial", 10, "bold")).pack(pady=10)
        self.lbl_task = tk.Label(self.top, text="Initializing...", font=("Arial", 9))
        self.lbl_task.pack(pady=5)

        self.progress = ttk.Progressbar(self.top, mode='indeterminate')
        self.progress.pack(fill='x', padx=20, pady=10)
        self.progress.start(10)

    def update_text(self, text):
        self.lbl_task.config(text=text)

    def close(self):
        self.top.destroy()

# --- LOGGING SYSTEM ---
LOG_FILE = "metadoon_session.log"

def safe_log(text):
    """
    Writes text to the GUI terminal and appends it to the log file.
    """
    # 1. Update GUI (Thread-safe)
    if 'terminal_output' in globals():
        terminal_output.after(0, lambda: _insert_log(text))
    
    # 2. Write to File
    try:
        # Open in append mode ('a'), create if doesn't exist
        with open(LOG_FILE, "a", encoding="utf-8") as f:
            # Add timestamp for file log only
            timestamp = datetime.now().strftime("[%Y-%m-%d %H:%M:%S] ")
            # Clean text (remove extra newlines for log file if needed)
            clean_text = text.strip()
            if clean_text:
                f.write(f"{timestamp}{clean_text}\n")
    except Exception as e:
        print(f"Failed to write to log file: {e}")

def _insert_log(text):
    terminal_output.config(state='normal')
    terminal_output.insert(tk.END, text)
    terminal_output.see(tk.END)
    terminal_output.config(state='disabled')

def _insert_log(text):
    terminal_output.config(state='normal')
    terminal_output.insert(tk.END, text)
    terminal_output.see(tk.END)
    terminal_output.config(state='disabled')

# --- CONFIG MANAGER ---
class ConfigContainer:
    def __init__(self):
        # REMOVED: stat_test and dist_method (Fixed in R)
        self.params = {
            "threads": 1,
            "fastq_maxdiffs": 30,
            "fastq_maxee": 1.0,
            "minuniquesize": 2,
            "clustering_method": "otu", 
            "id": 0.97,
            "sintax_cutoff": 0.8,
            "strand": "both",
            "use_default_chimera_db": True,
            "chimera_db": "",
            "16s_db_option": "rdp",
            "custom_16s_db": "",
            "color_palette": "viridis",
            "rarefaction_step": 100,
            "rarefaction_cex": 0.6,
            "enable_rarefaction": False,
            "rarefaction_depth": 1000,
            "abundance_top_n": 15,
            "core_top_n": 30
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
                    loaded = json.load(json_file)
                    self.params.update(loaded)
                return self.params
            else:
                return self.params
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load parameters: {str(e)}")
            return self.params

config_container = ConfigContainer()


def organize_output():
    """
    Organizes the messy Output folder into structured subfolders.
    """
    safe_log("\n>>> Organizing Output Files...\n")
    
    base_output = os.path.join(os.getcwd(), 'Output')
    if not os.path.exists(base_output):
        return

    # Define structure
    structure = {
        'Plots/Alpha_Diversity': ['alpha_diversity_', 'rarefaction_curve'],
        'Plots/Beta_Diversity': ['beta_diversity_'],
        'Plots/Taxonomy': ['abundance_', 'core_microbiome_'],
        'Plots/Differential_Abundance': ['deseq2_diff_', 'ancombc_diff_'],
        'Tables/Differential_Abundance': ['DESeq2_', 'ANCOMBC_'],
        'Tables/Statistics': ['permanova_result_']
    }

    # Create folders
    for folder in structure.keys():
        os.makedirs(os.path.join(base_output, folder), exist_ok=True)

    # Move files
    files = os.listdir(base_output)
    for f in files:
        src = os.path.join(base_output, f)
        if not os.path.isfile(src): continue
        
        moved = False
        for folder, patterns in structure.items():
            if any(p in f for p in patterns):
                dst = os.path.join(base_output, folder, f)
                shutil.move(src, dst)
                moved = True
                break
        
        # Optional: Move remaining files to a "Misc" folder if needed
        # if not moved: ...

    safe_log(">>> Output Organized Successfully.\n")
# --- CONFIG WINDOW ---
def configure_execution():
    config_window = tk.Toplevel()
    config_window.title("Configure Execution")
    config_window.geometry("450x650")
    config_window.resizable(True, True)

    main_frame = tk.Frame(config_window)
    main_frame.pack(fill="both", expand=True)

    canvas = tk.Canvas(main_frame)
    scrollbar = tk.Scrollbar(main_frame, orient="vertical", command=canvas.yview)
    scrollable_frame = tk.Frame(canvas)

    scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
    canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
    canvas.configure(yscrollcommand=scrollbar.set)

    # Scroll Handling
    def _on_mousewheel(event):
        try: canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        except tk.TclError: pass

    def _on_linux_scroll_up(event):
        try: canvas.yview_scroll(-1, "units")
        except tk.TclError: pass

    def _on_linux_scroll_down(event):
        try: canvas.yview_scroll(1, "units")
        except tk.TclError: pass
    
    config_window.bind_all("<MouseWheel>", _on_mousewheel)
    config_window.bind_all("<Button-4>", _on_linux_scroll_up)
    config_window.bind_all("<Button-5>", _on_linux_scroll_down)

    def on_close():
        config_window.unbind_all("<MouseWheel>")
        config_window.unbind_all("<Button-4>")
        config_window.unbind_all("<Button-5>")
        config_window.destroy()

    config_window.protocol("WM_DELETE_WINDOW", on_close)

    # Helper Functions
    def toggle_default_chimera_db():
        if default_chimera_db_var.get() == "yes":
            default_chimera_button.config(state=tk.DISABLED)
        else:
            default_chimera_button.config(state=tk.NORMAL)

    def upload_chimera_db():
        file_path = filedialog.askopenfilename(title="Select Chimera Database")
        if file_path: config_container.params["chimera_db"] = file_path

    def select_16s_db(event=None):
        option = db_16s_combobox.get()
        if option == "Use Custom 16S Database":
            file_path = filedialog.askopenfilename(title="Select 16S Database")
            if file_path:
                config_container.params["custom_16s_db"] = file_path
                config_container.params["16s_db_option"] = "custom"
        elif option == "Use RDP 16S Database":
            config_container.params["16s_db_option"] = "rdp"
        else:
            config_container.params["16s_db_option"] = "greengenes2"

    def save_parameters():
        try:
            cluster_map = {"OTU (Cluster 97%)": "otu", "ASV (Denoising/Unoising)": "asv"}
            
            params = {
                "threads": int(treads.get()),
                "fastq_maxdiffs": int(max_diffs_entry.get()),
                "fastq_maxee": float(max_ee_entry.get()),
                "minuniquesize": int(min_unique_size_entry.get()),
                "clustering_method": cluster_map.get(cluster_method_combobox.get(), "otu"),
                "id": float(id_entry.get()),
                "sintax_cutoff": float(sintax_cutoff_entry.get()),
                "strand": strand_option.get(),
                "use_default_chimera_db": default_chimera_db_var.get() == "yes",
                "chimera_db": config_container.params.get("chimera_db", ""),
                "16s_db_option": config_container.params.get("16s_db_option", "rdp"),
                "custom_16s_db": config_container.params.get("custom_16s_db", ""),
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
            messagebox.showerror("Input Error", f"Invalid input: {ve}")

    def load_metadata():
        files = filedialog.askopenfilenames(title="Select Metadata File (CSV/TSV)")
        if files:
            base_dir = os.getcwd()
            metadata_dir = os.path.join(base_dir, 'Metadata File')
            os.makedirs(metadata_dir, exist_ok=True)
            for f in os.listdir(metadata_dir): os.remove(os.path.join(metadata_dir, f))
            shutil.copy(files[0], metadata_dir)
            if 'terminal_output' in globals():
                try:
                    terminal_output.config(state='normal')
                    terminal_output.insert(tk.END, f"Metadata loaded: {os.path.basename(files[0])}\n")
                    terminal_output.config(state='disabled')
                except: pass

    def load_tree():
        files = filedialog.askopenfilenames(title="Select Tree File (.nwk)")
        if files:
            base_dir = os.getcwd()
            tree_dir = os.path.join(base_dir, 'Tree File')
            os.makedirs(tree_dir, exist_ok=True)
            for f in os.listdir(tree_dir): os.remove(os.path.join(tree_dir, f))
            dest_path = os.path.join(tree_dir, "tree.nwk")
            shutil.copy(files[0], dest_path)
            if 'terminal_output' in globals():
                try:
                    terminal_output.config(state='normal')
                    terminal_output.insert(tk.END, f"Tree loaded: {os.path.basename(files[0])}\n")
                    terminal_output.config(state='disabled')
                except: pass

    # --- GUI WIDGETS ---
    ttk.Label(scrollable_frame, text="--- VSEARCH Parameters ---", font=('Arial', 9, 'bold')).pack(pady=5)
    
    ttk.Label(scrollable_frame, text="Threads:").pack()
    treads = ttk.Combobox(scrollable_frame, values=list(range(1, ContainerGeneral.n_threads + 1)))
    treads.set(config_container.params.get("threads", 1))
    treads.pack(pady=2)

    ttk.Label(scrollable_frame, text="Max Diffs (Merge):").pack()
    max_diffs_entry = ttk.Entry(scrollable_frame)
    max_diffs_entry.insert(0, config_container.params.get("fastq_maxdiffs", 30))
    max_diffs_entry.pack(pady=2)

    ttk.Label(scrollable_frame, text="Max EE (Filter):").pack()
    max_ee_entry = ttk.Entry(scrollable_frame)
    max_ee_entry.insert(0, config_container.params.get("fastq_maxee", 1.0))
    max_ee_entry.pack(pady=2)

    ttk.Label(scrollable_frame, text="Min Unique Size:").pack()
    min_unique_size_entry = ttk.Entry(scrollable_frame)
    min_unique_size_entry.insert(0, config_container.params.get("minuniquesize", 2))
    min_unique_size_entry.pack(pady=2)

    ttk.Label(scrollable_frame, text="Analysis Type (Clustering):", font=('Arial', 9, 'bold')).pack(pady=(10,2))
    cluster_method_combobox = ttk.Combobox(scrollable_frame, values=["OTU (Cluster 97%)", "ASV (Denoising/Unoising)"], state="readonly")
    current_method = config_container.params.get("clustering_method", "otu")
    if current_method == "asv":
        cluster_method_combobox.set("ASV (Denoising/Unoising)")
    else:
        cluster_method_combobox.set("OTU (Cluster 97%)")
    cluster_method_combobox.pack(pady=2)

    ttk.Label(scrollable_frame, text="Identity % (Only for OTU):").pack()
    id_entry = ttk.Entry(scrollable_frame)
    id_entry.insert(0, config_container.params.get("id", 0.97))
    id_entry.pack(pady=2)

    ttk.Label(scrollable_frame, text="SINTAX Cutoff (0.0-1.0):").pack()
    sintax_cutoff_entry = ttk.Entry(scrollable_frame)
    sintax_cutoff_entry.insert(0, config_container.params.get("sintax_cutoff", 0.8))
    sintax_cutoff_entry.pack(pady=2)

    ttk.Label(scrollable_frame, text="Strand:").pack()
    strand_option = ttk.Combobox(scrollable_frame, values=["plus", "both"])
    strand_option.set(config_container.params.get("strand", "both"))
    strand_option.pack(pady=2)

    ttk.Separator(scrollable_frame, orient='horizontal').pack(fill='x', pady=10)
    default_chimera_db_var = tk.StringVar(value="yes" if config_container.params.get("use_default_chimera_db", True) else "no")
    tk.Checkbutton(scrollable_frame, text="Use Default Chimera DB (RDP Gold)", variable=default_chimera_db_var, onvalue="yes", offvalue="no", command=toggle_default_chimera_db).pack()
    default_chimera_button = tk.Button(scrollable_frame, text="Select Custom Chimera DB", command=upload_chimera_db, state=tk.DISABLED)
    default_chimera_button.pack(pady=2)

    ttk.Label(scrollable_frame, text="16S Database:").pack()
    db_16s_combobox = ttk.Combobox(scrollable_frame, values=["Use RDP 16S Database", "Use Greengenes2 16S Database", "Use Custom 16S Database"], state="readonly")
    current_db_opt = config_container.params.get("16s_db_option", "rdp")
    if current_db_opt == "greengenes2": db_16s_combobox.set("Use Greengenes2 16S Database")
    elif current_db_opt == "custom": db_16s_combobox.set("Use Custom 16S Database")
    else: db_16s_combobox.set("Use RDP 16S Database")
    db_16s_combobox.pack(pady=2)
    db_16s_combobox.bind("<<ComboboxSelected>>", select_16s_db)

    ttk.Separator(scrollable_frame, orient='horizontal').pack(fill='x', pady=10)
    ttk.Button(scrollable_frame, text="Upload Metadata File", command=load_metadata).pack(pady=5)
    ttk.Button(scrollable_frame, text="Upload Tree File (Optional)", command=load_tree).pack(pady=5)

    ttk.Separator(scrollable_frame, orient='horizontal').pack(fill='x', pady=10)
    ttk.Label(scrollable_frame, text="--- R Analysis Parameters ---", font=('Arial', 9, 'bold')).pack(pady=5)

    # REMOVED STAT TEST AND DIST METHOD SELECTORS

    ttk.Label(scrollable_frame, text="Color Palette:").pack()
    color_palette_combobox = ttk.Combobox(scrollable_frame, values=["viridis", "plasma", "inferno", "magma", "cividis", "wesanderson", "RColorBrewer"], state="readonly")
    color_palette_combobox.set(config_container.params.get("color_palette", "viridis"))
    color_palette_combobox.pack(pady=2)

    ttk.Label(scrollable_frame, text="Rarefaction Settings:", font=('Arial', 8, 'bold')).pack(pady=(10,2))
    enable_rarefaction_var = tk.BooleanVar(value=config_container.params.get("enable_rarefaction", False))
    tk.Checkbutton(scrollable_frame, text="Enable Rarefaction", variable=enable_rarefaction_var).pack()
    
    ttk.Label(scrollable_frame, text="Depth:").pack()
    rarefaction_depth_entry = ttk.Entry(scrollable_frame)
    rarefaction_depth_entry.insert(0, str(config_container.params.get("rarefaction_depth", 1000)))
    rarefaction_depth_entry.pack(pady=2)

    ttk.Label(scrollable_frame, text="Curve Step:").pack()
    rarefaction_step_entry = ttk.Entry(scrollable_frame)
    rarefaction_step_entry.insert(0, str(config_container.params.get("rarefaction_step", 100)))
    rarefaction_step_entry.pack(pady=2)

    ttk.Label(scrollable_frame, text="Curve Scale (cex):").pack()
    rarefaction_cex_entry = ttk.Entry(scrollable_frame)
    rarefaction_cex_entry.insert(0, str(config_container.params.get("rarefaction_cex", 0.6)))
    rarefaction_cex_entry.pack(pady=2)

    ttk.Label(scrollable_frame, text="Plot Limits (Top N):", font=('Arial', 8, 'bold')).pack(pady=(10,2))
    ttk.Label(scrollable_frame, text="Abundance Top N:").pack()
    abundance_top_n_entry = ttk.Entry(scrollable_frame)
    abundance_top_n_entry.insert(0, str(config_container.params.get("abundance_top_n", 15)))
    abundance_top_n_entry.pack(pady=2)

    ttk.Label(scrollable_frame, text="Core Top N:").pack()
    core_top_n_entry = ttk.Entry(scrollable_frame)
    core_top_n_entry.insert(0, str(config_container.params.get("core_top_n", 30)))
    core_top_n_entry.pack(pady=2)

    ttk.Button(scrollable_frame, text="Save Parameters", command=save_parameters, width=30).pack(pady=20)

    canvas.pack(side="left", fill="both", expand=True)
    scrollbar.pack(side="right", fill="y")

# --- FILE HANDLING ---
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

def generate_stats():
    terminal_output.config(state='normal')
    for file in container.files:
        if file not in container.stats:
            try:
                terminal_output.insert(tk.END, f"Calculating stats for {os.path.basename(file)}...\n")
                command = ["vsearch", "--fastq_eestats2", file, "--output", "-"]
                result = subprocess.run(command, check=True, capture_output=True, text=True)
                container.add_stats(file, result.stdout)
            except Exception as e:
                terminal_output.insert(tk.END, f"Error stats: {e}\n")
    terminal_output.config(state='disabled')

def display_stats():
    stats_text.config(state='normal')
    stats_text.delete(1.0, tk.END)
    for file in container.files:
        if file in container.stats:
            stats_text.insert(tk.END, f"--- {os.path.basename(file)} ---\n")
            stats_text.insert(tk.END, f"{container.stats[file]}\n\n")
    stats_text.config(state='disabled')

def open_file():
    files = filedialog.askopenfilenames(title="Select FASTQ Files",
                                        filetypes=(("FASTQ Files", "*.fastq"), ("All Files", "*.*")))
    for file in files:
        container.add_file(file)
        file_list.insert(tk.END, os.path.basename(file))

    generate_stats()
    display_stats()

# --- WORKSPACE CLEANUP ---
def clean_workspace():
    safe_log("\n>>> Starting Workspace Cleanup...\n")
    base_dir = os.getcwd()
    
    folders_to_remove = ['Merged', 'Filtered', 'Dereplicated', 'OTUs', 'Taxonomy', 'DB', 'FullFiles', 'Output', 'Metadata File', 'Tree File']
    files_to_remove = ['pipeline_params.json', 'Metadoon_Report.html', 'Rplots.pdf', '.RData', '.Rhistory','metadoon_session.log']

    for folder in folders_to_remove:
        dir_path = os.path.join(base_dir, folder)
        if os.path.exists(dir_path):
            try:
                shutil.rmtree(dir_path)
                safe_log(f"Deleted folder: {folder}\n")
            except Exception as e:
                safe_log(f"Error deleting folder {folder}: {e}\n")

    for file in files_to_remove:
        file_path = os.path.join(base_dir, file)
        if os.path.exists(file_path):
            try:
                os.remove(file_path)
                safe_log(f"Deleted file: {file}\n")
            except Exception as e:
                safe_log(f"Error deleting file {file}: {e}\n")

    try:
        container.files.clear()
        container.stats.clear()
        file_list.delete(0, tk.END)
        stats_text.config(state='normal')
        stats_text.delete(1.0, tk.END)
        stats_text.config(state='disabled')
        safe_log(">>> Workspace cleaned successfully. Ready for new analysis.\n")
    except Exception as e:
        safe_log(f"Error refreshing UI: {e}\n")

# --- PIPELINE EXECUTION ---
def run_pipeline():
    loading = LoadingWindow(root, "Running Pipeline")
    t = threading.Thread(target=execute_pipeline, args=(loading,))
    t.start()

def organize_output():
    """
    Organizes the messy Output folder into structured subfolders.
    """
    safe_log("\n>>> Organizing Output Files...\n")
    
    base_output = os.path.join(os.getcwd(), 'Output')
    if not os.path.exists(base_output):
        return

    # Define structure
    structure = {
        'Plots/Alpha_Diversity': ['alpha_diversity_', 'rarefaction_curve'],
        'Plots/Beta_Diversity': ['beta_diversity_'],
        'Plots/Taxonomy': ['abundance_', 'core_microbiome_'],
        'Plots/Differential_Abundance': ['deseq2_diff_', 'ancombc_diff_'],
        'Tables/Differential_Abundance': ['DESeq2_', 'ANCOMBC_'],
        'Tables/Statistics': ['permanova_result_']
    }

    # Create folders
    for folder in structure.keys():
        os.makedirs(os.path.join(base_output, folder), exist_ok=True)

    # Move files
    files = os.listdir(base_output)
    for f in files:
        src = os.path.join(base_output, f)
        if not os.path.isfile(src): continue
        
        moved = False
        for folder, patterns in structure.items():
            if any(p in f for p in patterns):
                dst = os.path.join(base_output, folder, f)
                shutil.move(src, dst)
                moved = True
                break
        
        # Optional: Move remaining files to a "Misc" folder if needed
        # if not moved: ...

    safe_log(">>> Output Organized Successfully.\n")

def execute_pipeline(loading_window):
    try:
        # Clear GUI Terminal
        terminal_output.after(0, lambda: terminal_output.delete(1.0, tk.END))
        
        # --- NEW: Clear Log File for new run ---
        if os.path.exists(LOG_FILE):
            os.remove(LOG_FILE)
            
        safe_log(">>> Initializing Metadoon Pipeline...\n")
        time.sleep(1.0)

        base_dir = os.getcwd()
        dirs = ['Merged', 'Filtered', 'Dereplicated', 'OTUs', 'Taxonomy', 'DB', 'FullFiles', 'Output', 'TempInput']
        for d in dirs:
            os.makedirs(os.path.join(base_dir, d), exist_ok=True)

        params = config_container.load_params()
        if not params:
            safe_log("Error: Could not load parameters.\n")
            return

        # ------------------ DOWNLOAD DBS ------------------
        loading_window.update_text("Checking Databases...")
        time.sleep(0.5)
        
        # Chimera
        if params["use_default_chimera_db"]:
            chimera_db = os.path.join(base_dir, 'DB', 'rdp_gold.fa')
            if not os.path.exists(chimera_db):
                safe_log("Downloading Chimera DB...\n")
                subprocess.run(["wget", "http://drive5.com/uchime/rdp_gold.fa", "-O", chimera_db, "-q"], check=True)
        else:
            chimera_db = params["chimera_db"]

        # SINTAX 16S
        if params["16s_db_option"] == "rdp":
            sintax_db = os.path.join(base_dir, 'DB', 'rdp_16s_v16.fa.gz')
            if not os.path.exists(sintax_db):
                safe_log("Downloading RDP DB...\n")
                subprocess.run(["wget", "https://www.drive5.com/sintax/rdp_16s_v16.fa.gz", "-O", sintax_db, "-q"], check=True)
        elif params["16s_db_option"] == "greengenes2":
            sintax_db = os.path.join(base_dir, 'DB', '2024.09.seqs.fna.gz')
            if not os.path.exists(sintax_db):
                safe_log("Downloading Greengenes2 DB...\n")
                subprocess.run(["wget", "http://ftp.microbio.me/greengenes_release/current/2024.09.seqs.fna.gz", "-O", sintax_db, "-q"], check=True)
        else:
            sintax_db = params["custom_16s_db"]

        # ------------------ MERGE PAIRS ------------------
        loading_window.update_text("Merging Paired Reads...")
        merged_files = []
        
        for file in container.files:
            if "_R1_" in file:
                temp_input_dir = os.path.join(base_dir, 'TempInput')
                
                file_name = os.path.splitext(os.path.basename(file))[0].replace("_R1_", "_")
                reverse_file_orig = file.replace('_R1_', '_R2_')
                local_r1 = os.path.join(temp_input_dir, os.path.basename(file))
                local_r2 = os.path.join(temp_input_dir, os.path.basename(reverse_file_orig))
                
                shutil.copy2(file, local_r1)
                if os.path.exists(reverse_file_orig):
                    shutil.copy2(reverse_file_orig, local_r2)
                else:
                    safe_log(f"Error: Reverse file not found for {file_name}\n")
                    continue

                merged_output = os.path.join(base_dir, 'Merged', f"{file_name}_merged.fastq")
                
                cmd = [
                    "vsearch", "--fastq_mergepairs", local_r1, "--reverse", local_r2,
                    "--fastq_maxdiffs", str(params['fastq_maxdiffs']),
                    "--fastqout", merged_output, "--relabel", f"{file_name}_seq_"
                ]
                
                safe_log(f"Merging {file_name}...\n")
                proc = subprocess.run(cmd, capture_output=True, text=True)
                
                if proc.returncode != 0:
                    safe_log(f"VSEARCH Error on {file_name}:\n{proc.stderr}\n")
                    raise Exception("VSEARCH Merge failed. See log.")
                
                merged_files.append(merged_output)

        if not merged_files:
            safe_log("Error: No merged files created. Check input filenames (must contain _R1_).\n")
            return

        # Concatenate
        all_merged = os.path.join(base_dir, 'FullFiles', 'all_samples.fastq')
        with open(all_merged, 'wb') as outfile:
            for fname in merged_files:
                with open(fname, 'rb') as infile:
                    shutil.copyfileobj(infile, outfile)

        # ------------------ FILTERING ------------------
        loading_window.update_text("Filtering Reads...")
        safe_log("Filtering reads...\n")
        time.sleep(0.5)
        all_filtered = os.path.join(base_dir, 'Filtered', 'all_filtered.fasta')
        cmd_filter = [
            "vsearch", "--fastq_filter", all_merged,
            "--fastq_maxee", str(params['fastq_maxee']),
            "--fastaout", all_filtered, "--fasta_width", "0"
        ]
        proc = subprocess.run(cmd_filter, capture_output=True, text=True)
        if proc.returncode != 0: raise Exception(f"Filtering failed: {proc.stderr}")

        # ------------------ DEREPLICATION ------------------
        loading_window.update_text("Dereplicating Sequences...")
        safe_log("Dereplicating...\n")
        time.sleep(0.5)
        all_derep = os.path.join(base_dir, 'Dereplicated', 'all_derep.fasta')
        proc = subprocess.run(["vsearch", "--derep_fulllength", all_filtered, "--output", all_derep, "--sizeout"], capture_output=True, text=True)
        if proc.returncode != 0: raise Exception(f"Dereplication failed: {proc.stderr}")

        # ------------------ CLUSTERING / DENOISING ------------------
        centroids = os.path.join(base_dir, 'OTUs', 'centroids.fasta')
        method = params.get("clustering_method", "otu")

        if method == "asv":
            loading_window.update_text("Denoising (Generating ASVs)...")
            safe_log("Denoising sequences (Unoising) to generate ASVs...\n")
            time.sleep(0.5)
            
            cmd = [
                "vsearch", "--cluster_unoise", all_derep,
                "--centroids", centroids,
                "--minsize", str(params['minuniquesize']),
                "--sizein", "--sizeout",
                "--relabel", "ASV_"
            ]
        else:
            loading_window.update_text("Clustering OTUs (97%)...")
            safe_log("Clustering OTUs (97% Identity)...\n")
            time.sleep(0.5)
            
            cmd = [
                "vsearch", "--cluster_size", all_derep, 
                "--id", str(params['id']), 
                "--centroids", centroids, 
                "--sizein", "--sizeout", 
                "--relabel", "OTU_"
            ]
            
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0: raise Exception(f"Clustering/Denoising failed: {proc.stderr}")

        # ------------------ CHIMERA REMOVAL ------------------
        loading_window.update_text("Removing Chimeras...")
        safe_log("Removing Chimeras (De Novo)...\n")
        time.sleep(0.5)
        denovo_nonchim = os.path.join(base_dir, 'OTUs', 'denovo_nonchim.fasta')
        subprocess.run(["vsearch", "--uchime_denovo", centroids, "--nonchimeras", denovo_nonchim], check=True)

        safe_log("Removing Chimeras (Reference)...\n")
        final_otus = os.path.join(base_dir, 'OTUs', 'otus.fasta')
        proc = subprocess.run(["vsearch", "--uchime_ref", denovo_nonchim, "--db", chimera_db, "--nonchimeras", final_otus], capture_output=True, text=True)
        if proc.returncode != 0: raise Exception(f"Chimera ref failed: {proc.stderr}")

        # ------------------ OTU TABLE ------------------
        loading_window.update_text("Generating Abundance Table...")
        safe_log("Generating Abundance Table...\n")
        time.sleep(0.5)
        otutab = os.path.join(base_dir, 'OTUs', 'otutab.txt')
        proc = subprocess.run(["vsearch", "--usearch_global", all_filtered, "--db", final_otus, "--id", str(params['id']), "--otutabout", otutab], capture_output=True, text=True)
        if proc.returncode != 0: raise Exception(f"OTU Table generation failed: {proc.stderr}")

        # ------------------ TAXONOMY ------------------
        loading_window.update_text("Assigning Taxonomy...")
        safe_log("Assigning Taxonomy (SINTAX)...\n")
        time.sleep(0.5)
        tax_raw = os.path.join(base_dir, 'Taxonomy', 'taxonomy_raw.txt')
        tax_final = os.path.join(base_dir, 'Taxonomy', 'taxonomy.txt')
        
        proc = subprocess.run([
            "vsearch", "--sintax", final_otus, "--db", sintax_db,
            "--tabbedout", tax_raw, "--sintax_cutoff", str(params['sintax_cutoff']),
            "--strand", params['strand']
        ], capture_output=True, text=True)
        if proc.returncode != 0: raise Exception(f"SINTAX failed: {proc.stderr}")

        # --- CLEAN TAXONOMY FILE ---
        safe_log("Formatting Taxonomy Table...\n")
        header = "OTU ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n"
        
        with open(tax_raw, 'r') as infile, open(tax_final, 'w') as outfile:
            outfile.write(header)
            for line in infile:
                parts = line.strip().split('\t')
                otu_id = parts[0].split(';')[0]
                
                if len(parts) >= 2:
                    tax_str = parts[1]
                    clean_tax = re.sub(r"\([\d\.]+\)", "", tax_str)
                    levels = clean_tax.split(',')
                    clean_levels = [re.sub(r"^[a-z]:", "", lvl) for lvl in levels]
                    while len(clean_levels) < 7: clean_levels.append("NA")
                    clean_levels = clean_levels[:7]
                    outfile.write(f"{otu_id}\t{'\t'.join(clean_levels)}\n")
                else:
                    outfile.write(f"{otu_id}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n")

        safe_log(">>> VSEARCH Pipeline Finished.\n")
        
        # Clean Temp Files
        if os.path.exists(os.path.join(base_dir, 'TempInput')):
            shutil.rmtree(os.path.join(base_dir, 'TempInput'))

        # ------------------ RUN R ANALYSIS ------------------
        if shutil.which("Rscript"):
            loading_window.update_text("Running R Analysis...")
            safe_log(">>> Starting R Analysis...\n")
            time.sleep(1.0)
            
            process = subprocess.run(["Rscript", "Analise.R"], capture_output=True, text=True)
            if process.returncode == 0:
                safe_log(process.stdout)
                safe_log("\n>>> R Analysis Completed Successfully!\n")
                root.after(0, lambda: messagebox.showinfo("Success", "Pipeline Finished!"))
            else:
                safe_log("Error in R Analysis:\n")
                safe_log(process.stderr)
                root.after(0, lambda: messagebox.showerror("R Error", "Analysis failed."))
        else:
            safe_log("Rscript not found.\n")

    except Exception as e:
        safe_log(f"\nCRITICAL ERROR: {str(e)}\n")
        print(traceback.format_exc())
    finally:
        organize_output()
        root.after(0, loading_window.close)


def Generate_report():
    loading = LoadingWindow(root, "Generating Report")
    t = threading.Thread(target=execute_report, args=(loading,))
    t.start()

def execute_report(loading_window):
    try:
        loading_window.update_text("Compiling HTML Report...")
        safe_log(">>> Generating HTML Report...\n")
        time.sleep(1.0)
        
        process = subprocess.run(["Rscript", "generate_report.R"], capture_output=True, text=True)
        
        if process.returncode == 0:
            safe_log("Report generated!\n")
            report_path = os.path.join(os.getcwd(), "Metadoon_Report.html")
            if os.path.exists(report_path):
                import webbrowser
                webbrowser.open(f"file://{report_path}")
        else:
            safe_log(f"Error generating report:\n{process.stderr}\n")
            
    except Exception as e:
        safe_log(f"Error: {str(e)}\n")
    finally:
        root.after(0, loading_window.close)

def save_analysis_results():

    parent_dir = filedialog.askdirectory(title="Select Parent Folder to Save Results")
    if not parent_dir:
        return

    try:

        timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        folder_name = f"Metadoon_Results_{timestamp}"
        dest_path = os.path.join(parent_dir, folder_name)
        
        os.makedirs(dest_path, exist_ok=True)

        terminal_output.config(state='normal')
        terminal_output.insert(tk.END, f"\n>>> Creating result folder: {folder_name}\n")
        terminal_output.config(state='disabled')
        
        curr = os.getcwd()
        

        targets = [
            'Output',             # Gráficos e Tabelas Estatísticas
            'Metadoon_Report.html', # Relatório Final
            'pipeline_params.json', # Parâmetros usados
            'OTUs',               # Tabelas de OTU e Sequências Finais
            'Taxonomy',           # Classificação Taxonômica
            'Metadata File',      # IMPORTANTE: Salva o metadado usado
            'Tree File',          # IMPORTANTE: Salva a árvore usada
            'Merged',             # Opcional: Reads unidos
            'Filtered',           # Opcional: Reads filtrados
            'metadoon_session.log'
        ]
        
        saved_any = False
        
        for item in targets:
            src = os.path.join(curr, item)
            if os.path.exists(src):
                final_dst = os.path.join(dest_path, item)
                
                
                if os.path.isdir(src):
                    if os.path.exists(final_dst): shutil.rmtree(final_dst)
                    shutil.copytree(src, final_dst)
               
                else:
                    shutil.copy2(src, final_dst)
                
                saved_any = True
        
        if saved_any:
            
            if messagebox.askyesno("Cleanup", f"Results saved in:\n{folder_name}\n\nDo you want to clean the workspace (delete temporary files)?"):
                clean_workspace()
            else:
                messagebox.showinfo("Saved", f"Results saved successfully in:\n{dest_path}")
                
                
                try:
                    if platform.system() == "Windows":
                        os.startfile(dest_path)
                    elif platform.system() == "Darwin": # macOS
                        subprocess.Popen(["open", dest_path])
                    else: # Linux
                        subprocess.Popen(["xdg-open", dest_path])
                except:
                    pass
        else:
            messagebox.showwarning("Warning", "No result files were found to save. Did you run the pipeline?")
            try: os.rmdir(dest_path) 
            except: pass

    except Exception as e:
        messagebox.showerror("Error", f"Failed to save: {e}")
        print(traceback.format_exc())

# --- MAIN ---
root = tk.Tk()
root.title("Metadoon - Microbiome Pipeline")
root.geometry("900x650")
system_os = platform.system()
if system_os == "Windows":
    root.iconbitmap("Metadoon.ico")
elif system_os == "Darwin":
    img = ImageTk.PhotoImage(file="Metadoon.icns")
    root.iconphoto(True, img)
else:  
    img = ImageTk.PhotoImage(file="Metadoon.png")
    root.iconphoto(True, img)
check_dependencies()

# Layout
left_frame = tk.Frame(root, width=300, bg="#f0f0f0")
left_frame.pack(side="left", fill="y")
right_frame = tk.Frame(root, bg="white")
right_frame.pack(side="right", fill="both", expand=True)

# Controls
tk.Label(left_frame, text="METADOON", font=("Arial", 16, "bold"), bg="#f0f0f0").pack(pady=20)
tk.Button(left_frame, text="1. Load FASTQ Files", command=open_file, width=25, height=2).pack(pady=5)
tk.Button(left_frame, text="2. Configure Parameters", command=configure_execution, width=25, height=2).pack(pady=5)
tk.Button(left_frame, text="3. RUN PIPELINE", command=run_pipeline, width=25, height=2, bg="lightblue").pack(pady=20)
tk.Button(left_frame, text="4. Generate Report", command=Generate_report, width=25, height=2).pack(pady=5)
tk.Button(left_frame, text="5. Save Results", command=save_analysis_results, width=25, height=2).pack(pady=5)

# Outputs
tk.Label(right_frame, text="Loaded Files", bg="white", font=("Arial", 10, "bold")).pack(anchor="w", padx=10, pady=(10,0))
file_list = tk.Listbox(right_frame, height=8)
file_list.pack(fill="x", padx=10, pady=5)

tk.Label(right_frame, text="File Statistics", bg="white", font=("Arial", 10, "bold")).pack(anchor="w", padx=10)
stats_text = tk.Text(right_frame, height=8, state='disabled')
stats_text.pack(fill="x", padx=10, pady=5)

tk.Label(right_frame, text="Pipeline Log", bg="white", font=("Arial", 10, "bold")).pack(anchor="w", padx=10)
terminal_output = tk.Text(right_frame, height=12, bg="#2b2b2b", fg="#00ff00")
terminal_output.pack(fill="both", expand=True, padx=10, pady=10)

root.mainloop()
