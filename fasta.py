import os 
import matplotlib.pyplot as plt  
from Bio import SeqIO 
from rich import print
import math

print("Welcome to Fread!")
file_path = input("Write the road to your file FASTA/FASTQ: ").strip()
                                                                                   
def read_sequences(file_path):
    if file_path.endswith(".fasta") or file_path.endswith(".fa"):               #Reading biological data
        file_format = "fasta" 
    elif file_path.endswith(".fastq"):
        file_format = "fastq"
    else:
        raise ValueError("File format is not supported!")

    return list(SeqIO.parse(file_path, file_format))
print()

sequences = []
try: 
    sequences = read_sequences(file_path)
    print(f"[bold]{len(sequences)} sequences were read successfully:[/bold]")   
    for record in sequences:                                                 
        print(f"{record.id}: {len(record.seq)} bp")
except Exception as e:
    print("Error:", e)
    exit()                                                                      
print()

if sequences:                                                                   # stats only if file is not empty
    print("[bold]                    ***** STATISTICAL ANALYSIS *****[/bold]")

    print(f"Nomber of sequences: {len(sequences)}")
    total_lenth = 0

    for record in sequences:
        total_lenth += len(record.seq)
    avarage_lenth = total_lenth / len(sequences)

    print(f"Average lenth = {avarage_lenth:.2f}bp")

    print("Contents of (GC):")

    total_gc = 0
    total_bases = 0 
    for record in sequences:
        seq = record.seq.upper()
        g_count = seq.count("G")
        c_count = seq.count("C")
        gc_count = g_count + c_count
        total_gc += gc_count
        total_bases += len(seq)
        gc_content = (gc_count / len(seq)) * 100
        print(f"     {record.id} = {gc_content:.2f}%")

    total_gc_content = (total_gc / total_bases) * 100   
    print(f"     Total = {total_gc_content:.2f}%")

    gc_count = g_count + c_count
    total_gc += gc_count
    total_bases += len(seq)
    gc_content = (gc_count / len(seq)) * 100
    print(f"     {record.id} = {gc_content:.2f}%")           #5 spaces at the beginning of the text

total_gc_content = (total_gc / total_bases) * 100   
print(f"     Total = {total_gc_content:.2f}%")

total_A = 0
total_T = 0
total_C = 0
total_G = 0 
total_bases = 0

for record in sequences:
    seq = record.seq.upper()
    a_count = seq.count("A")
    t_count = seq.count("T")
    c_count = seq.count("C")
    g_count = seq.count("G")

    total_A += a_count
    total_T += t_count
    total_C += c_count
    total_G += g_count
    total_bases += len(seq)

#FREQUENCIES in %
freq_A = (total_A / total_bases) * 100
freq_T = (total_T / total_bases) * 100
freq_C = (total_C / total_bases) * 100
freq_G = (total_G / total_bases) * 100
print()

print("[bold]                    ***** NUCLEOTIDE FREQUENCIES *****[/bold]")
print(f"A = {freq_A:2f}%")
print(f"T = {freq_T:2f}%")
print(f"C = {freq_C:2f}%")
print(f"G = {freq_G:2f}%")
print()

ids = [record.id for record in sequences]                              #VISUALISATION 
lengths = [len(record.seq) for record in sequences]

chunk_size = 5
num_chunks = math.ceil(len(sequences) / chunk_size)

fig, axes = plt.subplots(num_chunks, 1, figsize=(10, 5 * num_chunks))

if num_chunks == 1:
    axes = [axes]

for idx, ax in enumerate(axes):
    start = idx * chunk_size
    end = start + chunk_size
    chunk_ids = ids[start:end]
    chunk_lengths = lengths[start:end]

    if not chunk_ids:  
        ax.axis('off')
        continue

    bars = ax.bar(chunk_ids, chunk_lengths, color='skyblue', edgecolor='black')

    for bar, length in zip(bars, chunk_lengths):
        ax.text(bar.get_x() + bar.get_width() / 2, length + 5, str(length),
                 ha='center', va='bottom', fontsize=8)

    ax.set_title(f"Sequences {start+1}–{start+len(chunk_ids)}")
    ax.set_xlabel("Sequence ID")
    ax.set_ylabel("Length (bp)")

    max_len = max(chunk_lengths)
    step = 10 ** (len(str(max_len)) - 2)
    step = max(step, 200)
    ax.set_yticks(range(0, math.ceil(max_len / step) * step + step, step))

    ax.tick_params(axis='x', rotation=45)

fig.tight_layout()
plt.show()


