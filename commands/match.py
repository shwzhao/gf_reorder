#!/usr/bin/env python3

import sys

def read_paf(input_file):
    """读取 PAF 文件并返回 DataFrame"""
    import pandas as pd
    col_names = [
        "queryID", "queryLen", "queryStart", "queryEnd",
        "strand", "refID", "refLen", "refStart",
        "refEnd", "numResidueMatches", "lenAln", "mapQ"
    ]
    return pd.read_csv(input_file, sep='\t', header=None, names=col_names, usecols=range(12))


def filter_by_prefix(data, column, prefix):
    """根据指定列的前缀过滤数据"""
    if prefix:
        data = data[data[column].str.startswith(prefix)]
    return data


def filter_by_threshold(data, min_query_len, min_rev_comp_len):
    """根据 queryLenAgg 和 queryRevCompLen 的阈值过滤数据"""
    return data[
        (data["queryLenAgg"] >= min_query_len) &
        (data["queryRevCompLen"].abs() >= min_rev_comp_len)
    ]


def process_alignments(alignments, ref_prefix=None, query_prefix=None, min_query_len=0, min_rev_comp_len=0, uniq=False):
    """
    对齐数据进行处理并根据条件过滤：
    1. 根据前缀过滤 refID 和 queryID
    2. 计算 lenAlnStrand
    3. 按 queryID 和 refID 分组，计算 queryLenAgg 和 queryRevCompLen
    4. 选择 queryLenAgg 最大的记录
    5. 添加 queryRevComp 列
    6. 按 refID 进一步筛选
    7. 根据 queryLenAgg 和 queryRevCompLen 的阈值过滤
    """
    import pandas as pd
    import numpy as np

    # 根据前缀过滤 refID 和 queryID
    alignments = filter_by_prefix(alignments, "refID", ref_prefix)
    alignments = filter_by_prefix(alignments, "queryID", query_prefix)

    # 计算 lenAlnStrand
    alignments=alignments.copy()
    alignments.loc[:, "lenAlnStrand"] = np.where(alignments["strand"] == "-", -alignments["lenAln"], alignments["lenAln"])

    # 按 queryID 和 refID 分组
    grouped = alignments.groupby(["queryID", "refID"], as_index=False).agg(
        queryLenAgg=("lenAln", "sum"),
        queryRevCompLen=("lenAlnStrand", "sum")
    )

    # 选择 queryID 内 queryLenAgg 最大的行
    grouped = grouped.loc[grouped.groupby("queryID")["queryLenAgg"].idxmax()]

    # 添加 queryRevComp 列
    grouped["queryRevComp"] = np.where(grouped["queryRevCompLen"] > 0, "+", "-")

    if uniq:
        # 选择 refID 内 queryLenAgg 最大的行
        grouped = grouped.loc[grouped.groupby("refID")["queryLenAgg"].idxmax()]

    # 根据阈值过滤
    grouped = filter_by_threshold(grouped, min_query_len, min_rev_comp_len)

    return grouped


def save_to_output(data, output_file):
    """将处理后的数据保存为文件或输出到标准输出"""
    if output_file:
        data.to_csv(output_file, sep='\t', index=False)
    else:
        data.to_csv(sys.stdout, sep='\t', index=False)


def setup_parser(parser):
    match_parser = parser.add_parser('match', help='Process a PAF file to generate a summarized rename mapping TSV file.')
    match_parser.add_argument("-i", "--input", required=True, help="Input .paf filename")
    match_parser.add_argument("-o", "--output", help="Output .tsv filename (default: standard output)")
    match_parser.add_argument("-u", "--uniq", action="store_true", help="Get the best match of ref (default: False)")
    match_parser.add_argument("--ref-prefix", help="Filter by refID prefix (optional)")
    match_parser.add_argument("--query-prefix", help="Filter by queryID prefix (optional)")
    match_parser.add_argument("--min-query-len", type=int, default=0,
                        help="Minimum queryLenAgg value to retain (default: 0)")
    match_parser.add_argument("--min-rev-comp-len", type=int, default=0,
                        help="Minimum absolute queryRevCompLen value to retain (default: 0)")
    return match_parser


def run(args):
    # args = parse_arguments()
    alignments = read_paf(args.input)
    processed_data = process_alignments(
        alignments,
        ref_prefix=args.ref_prefix,
        query_prefix=args.query_prefix,
        min_query_len=args.min_query_len,
        min_rev_comp_len=args.min_rev_comp_len,
        uniq=args.uniq
    )
    save_to_output(processed_data, args.output)



def main():
    pass


if __name__ == "__main__":
    main()

