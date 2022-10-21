import argparse
from argparse import HelpFormatter
import re
from math import log

def set_purity(file: str) -> dict:
    """
    Set the tumor purity for each sample

    :param file: input file
    :type file: str
    :return: dict of sample purities
    :rtype: dict
    """
    purities = {}
    with open(file, 'r') as f:
        for line in f:
            sample, purity = line.strip().split('\t')
            purity = float(purity)
            if purity > 1:
                purity /= 100.0
            purities[sample] = purity
    return purities


def set_mad(file: str) -> dict:
    """
    Set the tumor MAD for each sample

    :param file: input file
    :return: dict of sample MADs
    """
    mads = {}
    with open(args.mad_file, 'r') as f:
        for line in f:
            sample, mad = line.strip().split('\t')
            mads[sample] = float(mad)
    return mads


def seq2c2fm(args) -> str:
    """
    Convert Seq2C log2gene output to OncoPrint output

    :param args: main program input argument
    :type args: argparse.Namespace
    :return: OncoPrint output
    :rtype: str
    """
    purity = {}
    if args.purity_file:
        purity = set_purity(args.purity_file)

    mad = {}
    if args.mad_file:
        mad = set_mad(args.mad_file)

    MINDEL = args.min_del
    MINAMP = args.min_amp  # 6 copies for pure sample
    MINEXONDEL = args.min_exon_del
    MINEXONAMP = args.min_exon_amp
    MINEXON = args.min_exon
    N = args.num  # If a breakpoint is called more than N samples, then it's deemed a false positive and filtered
    genes_gain = args.genes_gain.split(':')
    MAD = args.mad

    cols = ['Sample', 'Empty', 'Variant_Type', 'Gene', 'NA', '-', '-', 'Segment', '-', '-', 'Transf_LogRatio',
            'Segments', 'LogRatio', 'alteration', '-', '-', '-', '-', '-', '-', '-', 'Alteration']
    output = '\t'.join(cols) + '\n' if args.print_header else ''

    samples = {}
    with open(args.in_file, 'r') as f:
        for line in f:
            a = line.replace('\n', '').split('\t')
            sample, gene = a[0], a[1]
            if sample == 'Sample':  # skip header
                continue

            if args.reg:
                m = re.search(args.reg, sample)
                if m is not None:
                    sample = m.group(0)
                else:
                    continue
            samples[sample] = 1
            SMINDEL, SMINAMP, SMINEXONDEL, SMINEXONAMP = MINDEL, MINAMP, MINEXONDEL, MINEXONAMP

            pur = 1.0
            if sample in purity:
                pur = purity[sample]
            elif args.purity_percent:
                pur = args.purity_percent / 100.0 if args.purity_percent > 1 else args.purity_percent

            lr_str = a[6]
            lr = float(lr_str)

            # If MAD is available, CNV will only be called if the values exceed the background
            if sample in mad and abs(lr) <= MAD * mad[sample]:
                continue

            if args.purity_file or args.purity_percent:
                copy = (2.0 ** lr - 1 + pur) * 2.0 / pur;  # For tumor absolute cooy: lr = log2((N*p+2*(1-p))/2)
                if copy <= 0:  # to capture the cases where lr will be really small for homozygous deletions
                    lr = -10
                    lr_str = str(lr)
                else:
                    lr = log(copy/2)/log(2)

            desc = f"{a[10]} of {a[11]}" if a[8] == "BP" else f"{a[11]} of {a[11]}"
            if args.output_gain or gene in genes_gain:  # Only do whole gene for copy gains
                if lr >= 0.75 and lr < SMINAMP:
                    vals = [sample, "", "copy-number-alteration", a[1], "NA", "-", "-", f"{a[2]}:{a[3]}", "-", "-",
                            f"{2 ** lr * 2:.1f}", desc, lr_str, "gain", "-", "-", "-", "-", "-", "-", "-", "Gain"]
                    output += "\t".join(vals) + "\n"
                    continue

            _type = "Deletion" if lr < SMINDEL else "Amplification"

            if a[10] and a[8] == "BP" and (float(a[10]) >= MINEXON or (float(a[11]) - float(a[10])) >= MINEXON):
                lr_str = a[12]
                lr = float(lr_str)
                seg = ''
                m = re.match(r"(\S+)\(", a[14])
                if m is not None:
                    seg = m.group(1)
                    m2 = re.match(r"(\d+).*,(\d+)$", seg)
                    if m2 is not None:
                        seg = f"{m2.group(1)}-{m2.group(2)}"

                if a[9] == "Del":
                    _type = "Deletion"
                    desc = f"Del seg {seg}"
                elif a[9] == "Dup":
                    _type = "Duplication"
                    desc = f"Dup seg {seg}"
                    lr_str = a[13]
                    lr = float(lr_str)  # use the difference instead of absolute log2 ratio for duplications

                if not (lr >= SMINEXONAMP or lr <= SMINEXONDEL):
                    continue
            else:
                if not (lr >= SMINAMP or lr <= SMINDEL):
                    continue

            if a[15] and int(a[15]) >= N:
                continue

            if _type == "Duplication":
                vals = [sample, "", "rearrangement", a[1], "likely", "-", "-", f"{a[2]}:{a[3]}", "-", "-",
                        f"{2 ** lr * 2:.1f}", desc, lr_str, "-", a[1], a[1], desc, "-", "-", "-", "-", "Rearrangement"]
            else:
                vals = [sample, "", "copy-number-alteration", a[1], "NA", "-", "-", f"{a[2]}:{a[3]}", "-", "-",
                        f"{2 ** lr * 2:.1f}", desc, lr_str,
                        "loss" if _type == "Deletion" else "amplification", "-", "-", "-", "-", "-", "-", "-", _type]
            output += "\t".join(vals) + "\n"

    return output


def main():
    usage = '''seq2c2fm.py [-g] [-e exons] [-n reg] [-N num] [-A amp] [-a amp] [-D del] [-d del] [-p purity_file] [-P purity] -i lr2gene_output -o oncoPrint_output'''
    parser = argparse.ArgumentParser(prog='seq2c2fm', description='Convertor of Seq2C results to OncoPrint format',
                                     epilog=f'''This program will parse seq2c output and make calls for each gene and output in the format compatible with OncoPrint\nUsage: {usage}''',
                                     formatter_class=HelpFormatter)  # argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-k', '--print_header', action='store_true', help='Print header')
    parser.add_argument('-g', '--output_gain', action='store_true', default=False, help='Whether to output copy gains [4-5] copies')
    parser.add_argument('-p', '--purity_file', type=str,
                        help='A file that contains the tumor purity for all samples. Two columns, first is sample name, second is the purity in %% or fraction [0-1]')
    parser.add_argument('-P', '--purity_percent', type=float,
                        help='The purity. Default: 1 or 100%%, as is for cell lines. If set, all samples will assume to have the same purity')
    parser.add_argument('-n', '--reg', default=None,
                        help='The regular expression to extract sample names. Default: none')
    parser.add_argument('-m', '--mad_file', type=str,
                        help='A file contains the MAD values for all samples.  Two columns, first is the sample name, 2nd is the MAD in log2 fraction')
    parser.add_argument('-M', '--mad', type=float, default=3.0,
                        help='The MAD. Default: 3, as is for cell lines. If set, all samples will assume to have the same MAD')
    parser.add_argument('-N', '--num', type=int, default=5,
                        help='If an breakpoint has been called in >= num of samples, it is deemed false positive. Default: 5')
    parser.add_argument('-e', '--min_exon', type=int, default=1, help='Minimum number of exons/amplicon. Default: 1')
    parser.add_argument('-D', '--min_del', type=float, default=-2.0,
                        help='For whole gene: The log2ratio to determine that a gene is homozygously deleted. Default: -2.0')
    parser.add_argument('-A', '--min_amp', type=float, default=1.45,
                        help='For whole gene: The log2ratio to determine that a gene is amplified. Default: 1.45.')
    parser.add_argument('-d', '--min_exon_del', type=float, default=-2.0,
                        help='The minimum log2ratio to determine that 1-2 exons are deleted. Should be lower than [-d] to reduce false positives. Default: -2.5')
    parser.add_argument('-a', '--min_exon_amp', type=float, default=1.75,
                        help='The minimum log2ratio to determine that 1-2 exons are amplified. Should be larger than [-a] to reduce false positives. Default: 1.75')
    parser.add_argument('-G', '--genes_gain', type=str, default='MYC',
                        help='List of genes, seperated by ":", for which gain will be captured. Default: MYC')
    parser.add_argument('-i', '--in_file', required=True, help='input file that is the lr2gene output')
    parser.add_argument('-o', '--out_file', required=True, help='Output file')

    args = parser.parse_args()
    # print(args)
    output = seq2c2fm(args)
    with open(args.out_file, 'w') as fout:
        fout.write(output)


if __name__ == '__main__':
    main()
