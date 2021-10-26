#Copyright (c) 2018-2020 Gavin W. Wilson

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.


def snv2vcf_cmd(cargs):
    import numpy as NP
    import pandas as pd
    from .snvcmp import SampleData
    from .calc import write_MTX
    import sys, argparse, os
    import itertools

    parser = argparse.ArgumentParser(description='Create a minimal VCF and mtx matrices from SNV data by removing potential edits in REDIportal or overapping Alu elements')
    parser.add_argument('--max-homozygous', help='Percent of cells that only express the SNV must be less than this', type=float, default=99.9)
    parser.add_argument('--min-cells', help='Number of cells with coverage at a given SNV', type=int, default=10)
    parser.add_argument('--min-alts', help='Number of cells with a given SNV', type=int, default=5)
    parser.add_argument('--min-af', help='Minimum overall allele fraction (total alt / [total ref + alt])', type=float, default=0)
    parser.add_argument('-o', '--output', help='Output folder', type=str)

    parser.add_argument('-r', '--ref', help='Ref lengths text generated by scsnv index (to generate vcf header)', type=str)
    parser.add_argument('-f', '--fasta', help='Fasta reference file (for vcf header)', type=str)
    egroup = parser.add_mutually_exclusive_group(required=True)
    egroup.add_argument('-e', '--remove-edits', help='Filter strand corrected A-to-G changes that overlap REDIportal and/or Alu elements', action='store_true')
    egroup.add_argument('-k', '--keep-annotated', help='Only keep annotated edits (1000 genomes)', action='store_true')
    #egroup.add_argument('-n', '--only-edits', help='Filter strand corrected A-to-G changes that overlap REDIportal and/or Alu elements', action='store_true')
    #egroup.add_argument('-b', '--keep-unannotated', help='Only keep unannotated edits (1000 genomes)', action='store_true')

    group = parser.add_mutually_exclusive_group(required=True)
    #group.add_argument('-c', '--csv', help='Write CSV files compatible with scSplit', action='store_true')
    group.add_argument('-m', '--mtx', help='Write vartrix like MTX files compatible with souporcell / vireo', action='store_true')
    group.add_argument('-x', '--flamm', help='Write flammkuchen H5 file', action='store_true')

    parser.add_argument('annotated', help='Annotated mutation flammkuchen/deepdish file')
    args = parser.parse_args(cargs[1:])
    print(args)

    sdata = args.annotated
    gt_dir = args.output
    rdata = args.ref
    data = SampleData(sdata, 'scSNV')
    if not os.path.isdir(gt_dir):
        os.mkdir(gt_dir)

    tref = data.strand_ref
    talt = data.strand_alt


    snvs = data.snvs

    cell_counts = []
    alt_counts = []
    BB_counts = []
    AFs = []

    for s, e in zip(tref.indptr, tref.indptr[1:]):
        d1 = tref.data[s:e]
        d2 = talt.data[s:e]
        rsum = d1.sum()
        asum = d2.sum()
        AFs.append(100.0 * asum / (rsum + asum))
        alt_counts.append(NP.count_nonzero(d2))
        cell_counts.append(len(d2))
        BB_counts.append(100.0 * ((d2 > 0) & (d1 == 0)).sum() / len(d1)) 

    BB_counts = NP.array(BB_counts)
    cell_counts = NP.array(cell_counts)
    alt_counts = NP.array(alt_counts)
    AFs = NP.array(AFs)

    keep = NP.full(len(snvs), True, dtype='bool')

    if args.min_cells > 0:
        keep &= (cell_counts >= args.min_cells)
    if args.min_alts > 0:
        keep &= (alt_counts >= args.min_alts)
    if args.max_homozygous > 0:
        keep &= (BB_counts < args.max_homozygous)
    if args.min_af < 101:
        keep &= (AFs >= args.min_af)

    if args.remove_edits:
        snvs['repeat_type'] = ''
        snvs.loc[snvs.rep_family != '', 'repeat_type'] = 'Other'
        snvs.loc[snvs.rep_family.str.contains('Alu'), 'repeat_type'] = 'Alu'
        snvs['potential_edit'] = (snvs.strand_ref == 'A') & (snvs.strand_alt == 'G') & ((snvs.repeat_type == 'Alu') | snvs.REDIportal)
        if args.remove_edits:
            print("Removing edits")
            keep &= ~snvs['potential_edit']
    elif args.keep_annotated:
        print("Only keeping annotated SNVs")
        keep &= snvs.g1000

    fsnvs = snvs[keep].copy()

    print(f'Using {len(fsnvs)} / {len(snvs)} SNVs for genotyping')

    ridx = fsnvs.snv_idx.values.copy()
    tref = tref[ridx]
    talt = talt[ridx]
    if args.mtx:
        print(tref.data.dtype, tref.indptr.dtype, tref.indices.dtype, talt.data.dtype)
        print(__file__)
        write_MTX(tref, talt, os.path.join(gt_dir, f'refs.mtx'), os.path.join(gt_dir, f'alts.mtx'))
        print("Testing")
        #mmdebug(os.path.join(gt_dir, f'refs.mtx'))
        #mmdebug(os.path.join(gt_dir, f'alts.mtx'))
        if not os.path.isfile(f'{gt_dir}/barcodes.txt'):
            fout = open(f'{gt_dir}/barcodes.txt', 'w')
            for b in data.barcodes:
                fout.write(b + "\n")
            fout.close()
    elif args.csv:
        tref = tref.toarray()
        rows = [f'{i[0]}:{i[1] + 1}' for i in fsnvs.index]
        print("Writing reference allele csv")
        df = pd.DataFrame(data=tref, columns=data.barcodes, index=rows)
        df.to_csv(f'{gt_dir}/refs.csv')

        talt = talt.toarray()
        print("Writing alternative allele csv")
        df = pd.DataFrame(data=talt, columns=data.barcodes, index=rows)
        df.to_csv(f'{gt_dir}/alts.csv')

    elif args.flamm:
        print("not supporting yet")

    rdata = pd.read_csv(rdata, sep='\t', names=['chrom', 'rlen', 'comment'])
    print("Writing VCF")
    with open(f'{gt_dir}/snvs.vcf', 'w') as fp:
        fp.write('##fileformat=VCFv4.2\n')
        fp.write(f'##reference={args.fasta}\n')
        for r in rdata.itertuples():
            fp.write(f'##contig=<ID={r.chrom},length={r.rlen}>\n')
        fp.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for i, r in fsnvs.iterrows():
            fp.write(f'{i[0]}\t{i[1] + 1}\tscsnv_idx_{r.snv_idx}\t{r.ref}\t{r.alt}\t.\tPASS\t.\n')

    print("Writing SNV csv file")
    fsnvs.to_csv(f'{gt_dir}/snvs.csv')
