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

    # Frame principal que conterá a barra de rolagem e o conteúdo
    main_frame = tk.Frame(config_window)
    main_frame.pack(fill="both", expand=True)

    # Canvas para permitir a rolagem
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

    # Funções internas
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


    # Configurações de interface gráfica dentro do frame rolável
    ttk.Label(scrollable_frame, text="Threads:").pack(pady=5)
    opcoes = list(range(1, ContainerGeneral.n_threads + 1))
    treads = ttk.Combobox(scrollable_frame)
    treads['values'] = opcoes
    treads.set(opcoes[0])
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

    # Opções de banco de dados Chimera
    default_chimera_db_var = tk.StringVar(value="yes")
    tk.Checkbutton(scrollable_frame, text="Use Default Chimera Database", variable=default_chimera_db_var, onvalue="yes", offvalue="no", command=toggle_default_chimera_db).pack(pady=5)
    default_chimera_button = tk.Button(scrollable_frame, text="Select Chimera Database", command=upload_chimera_db, state=tk.DISABLED)
    default_chimera_button.pack(pady=5)

    # Opções de banco de dados 16S
    ttk.Label(scrollable_frame, text="16S Database Options:").pack(pady=5)
    db_16s_combobox = ttk.Combobox(scrollable_frame, values=["Use RDP 16S Database", "Use Greengenes2 16S Database", "Use Custom 16S Database"])
    db_16s_combobox.set("Use RDP 16S Database")
    db_16s_combobox.pack(pady=5)
    db_16s_combobox.bind("<<ComboboxSelected>>", lambda e: select_16s_db())
    #Upload metadata
    ttk.Label(scrollable_frame, text="Metadata File:").pack(pady=5)
    ttk.Label(scrollable_frame, text="Metadata file here", background="lightgreen", foreground="red").pack(padx=2,pady=2)
    metadata_entry = ttk.Button(scrollable_frame, text="Upload Metadata File", command=load_metadata)
    metadata_entry.pack(pady=30)
    #Upload tree
    ttk.Label(scrollable_frame, text="Tree (Optional):").pack(pady=5)
    ttk.Label(scrollable_frame, text="Tree here", background="lightgreen", foreground="red").pack(padx=2,pady=2)
    tree_entry = ttk.Button(scrollable_frame, text="Upload Tree File", command=load_tree)
    tree_entry.pack(pady=30)
    # Botão para salvar os parâmetros configurados
    save_button = ttk.Button(scrollable_frame, text="Save Parameters", command=save_parameters)
    save_button.pack(pady=20)

    # Configurar e exibir o canvas e a barra de rolagem
    canvas.pack(side="left", fill="both", expand=True)
    scrollbar.pack(side="right", fill="y")

# Classe para armazenar arquivos carregados

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
                terminal_output.insert(tk.END, f"Gerando estatísticas para {os.path.basename(file)}...\n")
                command = ["vsearch", "--fastq_eestats2", file, "--output", "-"]
                result = subprocess.run(command, check=True, capture_output=True, text=True)
                container.add_stats(file, result.stdout)
                terminal_output.insert(tk.END, f"Estatísticas geradas para {os.path.basename(file)}\n")
            except FileNotFoundError:
                terminal_output.insert(tk.END, "vsearch não encontrado no sistema.\n")
                print("vsearch não encontrado no sistema.")
            except CalledProcessError as e:
                terminal_output.insert(tk.END, f"Erro ao calcular estatísticas para {file}:\n{e.stderr}\n")
                print(f"Erro ao calcular estatísticas para {file}:\n{e.stderr}\n")
            except Exception as e:
                terminal_output.insert(tk.END, f"Erro ao executar o vsearch: {e}\n{traceback.format_exc()}\n")
                print(f"Erro ao executar o vsearch: {e}\n{traceback.format_exc()}\n")
    terminal_output.config(state='disabled')

def display_stats():
    stats_text.delete(1.0, tk.END)
    for file in container.files:
        if file in container.stats:
            stats_text.insert(tk.END, f"Estatísticas de {os.path.basename(file)}:\n")
            stats_text.insert(tk.END, f"{container.stats[file]}\n")
        else:
            stats_text.insert(tk.END, f"Estatísticas para {file} não encontradas.\n")

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
            terminal_output.insert(tk.END, "Iniciando o Pipeline...\n")
        except:
            terminal_output.insert(tk.END, " .\n")

            # 1. Merge Pairs - Mesclando pares de leituras (R1 e R2)
            files_names = []
            merged_name_files = []
            derep_files = []
            clustered_OTUs = []
            # Diretórios de saída para os diferentes estágios do pipeline
            base_dir = os.getcwd()  # Use o diretório atual como base
            merged_dir = os.path.join(base_dir, 'Merged')
            filter_dir = os.path.join(base_dir, 'Filtered')
            derep_dir = os.path.join(base_dir, 'Dereplicated')
            otu_dir = os.path.join(base_dir, 'OTUs')
            tax_dir = os.path.join(base_dir, 'Taxonomy')
            db_dir = os.path.join(base_dir, 'DB')

            # Criação de diretórios conforme necessário
            os.makedirs(merged_dir, exist_ok=True)
            os.makedirs(filter_dir, exist_ok=True)
            os.makedirs(derep_dir, exist_ok=True)
            os.makedirs(otu_dir, exist_ok=True)
            os.makedirs(tax_dir, exist_ok=True)
            os.makedirs(db_dir, exist_ok=True)

            # Carregar parâmetros
            params = config_container.load_params()

            # Baixar banco de dados Chimera, se necessário
            if params["use_default_chimera_db"]:
                chimera_db = os.path.join(db_dir, 'rdp_gold.fa')
                if not os.path.exists(chimera_db):
                    terminal_output.insert(tk.END, "Downloading default Chimera DB...\n")
                    os.system(f"wget http://drive5.com/uchime/rdp_gold.fa --no-check-certificate -O {chimera_db}")
            else:
                chimera_db = params["chimera_db"]

            # Seleção do banco de dados 16S
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

           # 1. Mergin Reads Foward and Reverse
            for file in container.files:
                if "_R1_" in file:
                    file_name = os.path.splitext(os.path.basename(file))[0].replace("_R1_", "_")
                    reverse_file = file.replace('_R1_', '_R2_')
                    merged_output = os.path.join(merged_dir, f"{file_name}_merged.fastq")
                    cmd_merge = f"vsearch --fastq_mergepairs '{file}' --reverse '{reverse_file}' --fastq_maxdiffs 30 --fastq_maxdiffpct 10 --fastqout '{merged_output}' --relabel {file_name}_seq_ --fastq_eeout"
                    terminal_output.insert(tk.END, f"\nMerging paired reads for {file_name} and {reverse_file}...\n")
                    terminal_output.insert(tk.END, f"Runing:{cmd_merge}\n")
                    os.system(cmd_merge)
                    merged_name_files.append(merged_output)
                    files_names.append(file_name)
            for merged in merged_name_files:
                print("Merged Files: \n")
                print(merged)
             # Criar diretório para arquivos combinados de amostras
            os.makedirs("FullFiles", exist_ok=True)
            subprocess.run(f"cat {merged_dir}/*_merged.fastq > FullFiles/all_single_sample.fastq", shell=True, check=True)
            print("Arquivos combinados com sucesso.")
            #filtrando as leituras
            terminal_output.insert(tk.END, "\nFiltering reads...\n")
            cmd_filter = f"vsearch --fastq_filter 'FullFiles/all_single_sample.fastq' --fastq_maxee {params['fastq_maxee']} --sizein --sizeout --fastq_maxns 0 --fastaout {filter_dir}/all_filtered.fasta --fasta_width 0"
            terminal_output.insert(tk.END, f"Runing:\n\n {cmd_filter}\n")
            os.system(cmd_filter)
            terminal_output.insert(tk.END, "Reads filtered successfully.\n")

            # Mostrar contagem de sequências no nível da amostra
            print("\nSequences count at sample level:\n")
            os.system(f'cat {filter_dir}/all_filtered.fasta | grep -c "^>"')

            # Dereplicação adicional de todas as amostras combinadas
            cmd_derep_all = f"vsearch --derep_fulllength '{filter_dir}/all_filtered.fasta' --threads {params['threads']} --sizein --sizeout --output {derep_dir}/all_derep.fasta"
            os.system(cmd_derep_all)

            # Mostrar contagem de sequências desduplicadas
            print("Unique dereplicated sequences: \n")
            os.system(f"grep -c '^>' {derep_dir}/all_derep.fasta")

            # Clusterização para criar centroids
            cmd_cluster = f"vsearch --cluster_size '{derep_dir}/all_derep.fasta' --threads {params['threads']} --id {params['id']} --strand {params['strand']} --sizein --sizeout --centroids {otu_dir}/samples_centroids.fasta"
            os.system(cmd_cluster)

            # Ordenar os centroids e remover singletons
            cmd_sort = f"vsearch --sortbysize '{otu_dir}/samples_centroids.fasta' --threads {params['threads']} --sizein --sizeout --minsize {params['minuniquesize']} --output {otu_dir}/samples_centroids_sorted.fasta"
            os.system(cmd_sort)

            print("Non-singleton clusters: ")
            os.system(fr'grep -c "^>" {otu_dir}/samples_centroids_sorted.fasta')

            # Detecção de chimeras de novo
            cmd_uchime_denovo = f"vsearch --uchime_denovo '{otu_dir}/samples_centroids_sorted.fasta' --sizein --sizeout --fasta_width 0 --qmask none --nonchimeras {otu_dir}/denovo.nonchimeras.fasta"
            os.system(cmd_uchime_denovo)

            print("Unique sequences after de novo chimera detection:")
            os.system(fr'grep -c "^>" {otu_dir}/denovo.nonchimeras.fasta')

            # Detecção de chimeras baseada em referência
            chimera_db = params["chimera_db"] if not params["use_default_chimera_db"] else os.path.join(db_dir, "rdp_gold.fa")
            cmd_uchime_ref = f"vsearch --uchime_ref '{otu_dir}/denovo.nonchimeras.fasta' --threads {params['threads']} --db '{chimera_db}' --sizein --sizeout --fasta_width 0 --qmask none --dbmask none --nonchimeras {otu_dir}/nonchimeras.fasta"
            os.system(cmd_uchime_ref)

            print("Unique sequences after reference-based chimera detection:")
            os.system(f'grep -c "^>" {otu_dir}/nonchimeras.fasta')

            # Filtrar e relabelar OTUs
            cmd_filter_otus = f"vsearch --fastx_filter '{otu_dir}/nonchimeras.fasta' --threads {params['threads']} --sizein --sizeout --fasta_width 0 --relabel OTU_ --fastaout {otu_dir}/otus.fasta"
            os.system(cmd_filter_otus)

            # Criar tabela de OTUs
            cmd_otu_table = f"vsearch --usearch_global '{filter_dir}/all_filtered.fasta' --threads {params['threads']} --db '{otu_dir}/otus.fasta' --id {params['id']} --strand {params['strand']} --otutabout {otu_dir}/otutab.txt"
            os.system(cmd_otu_table)

            # Configurar base de dados SINTAX
            sintax_db = {
                "rdp": os.path.join(db_dir, "rdp_16s_v16.fa.gz"),
                "greengenes2": os.path.join(db_dir, "2024.09.seqs.fna.gz"),
                "custom": params["custom_16s_db"]
            }[params["16s_db_option"]]

            # Comando para análise de taxonomia usando SINTAX
            cmd_sintax = f"vsearch --sintax '{otu_dir}/otus.fasta' --db '{sintax_db}' --tabbedout {tax_dir}/taxonomy.txt --sintax_cutoff {params['sintax_cutoff']} --strand {params['strand']}"
            os.system(cmd_sintax)

            # Define o cabeçalho a ser adicionado
            header = "#OTU ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\t\n"

            # Caminho do diretório contendo o arquivo taxonomy.txt
            file_path = f"{tax_dir}/taxonomy.txt"

            # Verificar se o arquivo existe
            if os.path.exists(file_path):
                # Lista para armazenar o conteúdo processado
                processed_content = []

                # Ler e processar o conteúdo do arquivo original
                with open(file_path, 'r') as original_file:
                    for line in original_file:
                        # Divide a linha em colunas separadas por tabulação
                        columns = line.strip().split('\t')

                        # Remove tudo após o ponto e vírgula na primeira coluna
                        columns[0] = re.sub(r";.*", "", columns[0])

                        # Mantém apenas as três primeiras colunas
                        if len(columns) > 3:
                            columns = columns[:3]

                        # Substitui as vírgulas por tabulação e remove parênteses e seu conteúdo nas colunas após a primeira
                        columns = [columns[0]] + [re.sub(r"\(.*?\)", "", col.replace(',', '\t').replace('+', '').replace('d:', '').replace('p:', '').replace('c:', '').replace('o:', '').replace('f:', '').replace('g:', '').replace('s:', '')) for col in columns[1:]]

                        # Junta as colunas em uma linha modificada e adiciona à lista de conteúdo processado
                        processed_line = '\t'.join(columns) + '\n'
                        processed_content.append(processed_line)

                # Escrever o novo conteúdo com o cabeçalho adicionado
                with open(file_path, 'w') as modified_file:
                    modified_file.write(header)
                    modified_file.writelines(processed_content)

                print("Arquivo taxonomy.txt processado e cabeçalho adicionado com sucesso!")

            else:
                print(f"O arquivo {file_path} não foi encontrado.")

            try:
                subprocess.run(["Rscript", "Analise.R"], check=True, capture_output=True,
                               text=True)  # executa o script de analise
                terminal_output.insert(tk.END, "Análise concluída com sucesso.\n")
            except subprocess.CalledProcessError as e:
                terminal_output.insert(tk.END, f"Erro na execução da análise R:\n{e.stderr}\n")  # mostra erros detalhados
                return  # interrompe a execução caso ocorra um erro
            terminal_output.config(state='normal')
            terminal_output.insert(tk.END, "Pipeline Starting...\n")
    except Exception as e:
        terminal_output.insert(tk.END, f"Ocorreu um erro inesperado:\n")
        terminal_output.insert(tk.END, f"{traceback.format_exc()}\n")  # Imprime o traceback completo
        print(traceback.format_exc())  # Imprime o traceback no console também


#MAIN WINDOW
root = tk.Tk()
root.title("Metadoon Version 1.0")
# Verifique se é Windows ou não
if root.tk.call('tk', 'windowingsystem') == 'win32':
    root.iconbitmap("Metadoon.ico")
else:
    img = ImageTk.PhotoImage(file="./Metadoon.png")
    root.iconphoto(True, img)
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

# Dividindo a janela em duas colunas (Files e Estatísticas)
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
