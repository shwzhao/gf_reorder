#!/usr/bin/env python3
from Bio import SeqIO
from Bio.Seq import Seq

def reverse_coordinates(start, end, chrom_length):
    """反转基因坐标"""
    new_start = chrom_length - end + 1
    new_end = chrom_length - start + 1
    return new_start, new_end

def load_rename_mapping(rename_file, genome_file, prefix=None, length_threshold=None):
    """解析 rename 文件，构建名称映射表并应用筛选逻辑"""
    rename_dict = {}
    genome_lengths = {record.id: len(record.seq) for record in SeqIO.parse(genome_file, "fasta")}

    with open(rename_file, 'r') as f:
        for line in f:
            if line.startswith("queryID"):
                continue  # 跳过标题行
            parts = line.strip().split()
            queryID, refID, queryRevComp = parts[0], parts[1], parts[4]

            if queryID not in genome_lengths:
                continue
            chrom_length = genome_lengths[queryID]

            if length_threshold is not None:
                if chrom_length < length_threshold:
                    continue
            elif prefix is not None:
                if not queryID.startswith(prefix):
                    continue

            rename_dict[queryID] = {
                'new_name': refID,
                'length': chrom_length,
                'reverse': queryRevComp == '-'
            }
    return rename_dict

def update_gff(input_gff, rename_dict, output_gff, log_entries):
    """更新 GFF 文件"""
    with open(input_gff, 'r') as in_f, open(output_gff, 'w') as out_f:
        for line in in_f:
            line = line.strip()  # 去掉首尾的空白字符
            if not line:  # 跳过空行
                continue
            if line.startswith("#"):
                out_f.write(line + 'n')
                continue

            fields = line.strip().split('\t')
            if len(fields) < 8:  # 跳过字段不足的行
                # log_entries.append(f"Skipped invalid line: {line}")
                continue
            chrom, start, end, strand = fields[0], int(fields[3]), int(fields[4]), fields[6]

            if chrom in rename_dict:
                new_name = rename_dict[chrom]['new_name']
                chrom_length = rename_dict[chrom]['length']
                reverse = rename_dict[chrom]['reverse']

                if reverse:
                    start, end = reverse_coordinates(start, end, chrom_length)
                    strand = "-" if strand == "+" else "+"

                log_entries.append(f"GFF: {chrom} -> {new_name}, Reverse: {'Yes' if reverse else 'No'}, Length: {chrom_length}")
                chrom = new_name

            fields[0], fields[3], fields[4], fields[6] = chrom, str(start), str(end), strand
            out_f.write('\t'.join(fields) + '\n')

def update_fasta(genome_file, rename_dict, output_genome, log_entries):
    """更新基因组序列"""
    with open(output_genome, 'w') as out_f:
        for record in SeqIO.parse(genome_file, "fasta"):
            if record.id in rename_dict:
                new_name = rename_dict[record.id]['new_name']
                reverse = rename_dict[record.id]['reverse']
                length = rename_dict[record.id]['length']

                log_entries.append(f"FASTA: {record.id} -> {new_name}, Reverse: {'Yes' if reverse else 'No'}, Length: {length}")
                record.id = new_name
                record.description = new_name

                if reverse:
                    record.seq = record.seq.reverse_complement()
            SeqIO.write(record, out_f, "fasta")

def save_log(log_entries, log_file):
    """保存日志到文件或打印到屏幕"""
    if log_file:
        with open(log_file, 'w') as log_f:
            log_f.write("\n".join(log_entries) + "\n")
    else:
        print("\n".join(log_entries))

def setup_parser(parser):
    reorder_parser = parser.add_parser('reorder', help='Update GFF and FASTA files based on rename mapping.')
    reorder_parser.add_argument('-i', '--input_gff', required=True, help='Path to the input GFF file.')
    reorder_parser.add_argument('-g', '--input_genome', required=True, help='Path to the query genome FASTA file.')
    reorder_parser.add_argument('-r', '--rename_file', required=True, help='Path to the rename mapping file.')
    reorder_parser.add_argument('-I', '--output_gff', required=True, help='Path to the output GFF file.')
    reorder_parser.add_argument('-G', '--output_genome', required=True, help='Path to the output genome FASTA file.')
    reorder_parser.add_argument('--prefix', help='Prefix for queryID to include for renaming.')
    reorder_parser.add_argument('--length_threshold', type=int, help='Minimum length of sequences to include for renaming.')
    reorder_parser.add_argument('--log_file', help='Path to the output log file (default: print to screen).')
    return reorder_parser

def run(args):
    if args.length_threshold is not None:
        print(f"Using length_threshold: {args.length_threshold}")
    elif args.prefix is not None:
        print(f"Using prefix: {args.prefix}")
    else:
        print("Error: Either --length_threshold or --prefix must be provided.")
        return

    rename_dict = load_rename_mapping(args.rename_file, args.input_genome, prefix=args.prefix, length_threshold=args.length_threshold)

    log_entries = []
    update_gff(args.input_gff, rename_dict, args.output_gff, log_entries)
    update_fasta(args.input_genome, rename_dict, args.output_genome, log_entries)

    save_log(log_entries, args.log_file)

def main():
    pass

if __name__ == "__main__":
    main()

