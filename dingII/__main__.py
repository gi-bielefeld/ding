from argparse import ArgumentParser, FileType, ArgumentTypeError
import logging
import dingII.dingII_generate as dingII_generate
import sys
import dingII.dingII_parsesol as dingII_parsesol


def matching_range(v):
    iv = float(v)
    if iv < 0.0 or iv > 1.0:
        raise ArgumentTypeError("%s is an invalid range bound. Please only provide values in [0,1]." % v)
    return iv

def configure_logging():
    # Create logger
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    # Create console handler and set level to debug
    ch = logging.StreamHandler(sys.stderr)
    ch.setLevel(logging.DEBUG)

    # Create formatter
    formatter = logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s')

    # Add formatter to ch
    ch.setFormatter(formatter)

    # Add ch to logger
    logger.addHandler(ch)


def define_generate_subcommand(subparsers):
    generate_parser = subparsers.add_parser('generate',
                                            help='Generate the DING ILP or create a custom model file in order to fine tune how many genes per family are to be matched')
    group = generate_parser.add_mutually_exclusive_group()
    group.add_argument('-mm', '--maximal', action='store_true', help='Set matching model to maximal matching.')
    group.add_argument('-em', '--exemplary', action='store_true', help='Set matching model to exemplary matching.')
    group.add_argument('-im', '--intermediate', action='store_true',
                       help='Set matching model to intermediate matching.')
    group.add_argument('-r', '--range', type=matching_range, nargs=2,
                       help='Provide upper and lower percentiles to be matched per marker in range [0,1]. Actual discrete bounds will always be rounded up.')
    generate_parser.add_argument('-c', '--custom', type=FileType('r'), action='store',
                                 help='Provide a custom matching file.')
    add_unimog_parsing_groups(generate_parser)
    writewhat = generate_parser.add_mutually_exclusive_group(required=True)
    writewhat.add_argument('--writemodel', type=FileType('w'),
                           help='Write the matching model to a file in order to customize it.')
    writewhat.add_argument('--writeilp', type=FileType('w'), help='Write the resulting ILP to the specified file.')
    generate_parser.set_defaults(func=dingII_generate.main)


def define_parsesol_subcommand(subparsers):
    parsesol_parser = subparsers.add_parser('parsesol',
                                            help='Parse a gurobi solution into a distance and optionally give a matching')
    add_unimog_parsing_groups(parsesol_parser)
    parsesol_parser.add_argument('-m', '--matching', type=FileType('w'),
                                 help='Give the matching as a pair of indexed genomes.')
    g = parsesol_parser.add_mutually_exclusive_group(required=True)
    g.add_argument('--solgur', type=FileType('r'), help='Gurobi solution file with a single solution.')
    parsesol_parser.add_argument('--runs', type=FileType('w'),
                                 help='Write runs of indels to the specified file. Format: Each line represents a cycle, each consecutive sequence of oriented markers is a run. Runs in the same cycle are separated by tab characters an begin with an A-run or a tab character if no A-run exists.')
    parsesol_parser.add_argument('--numindels', action='store_true',
                                 help='Give a possible number of indels in the sorting scenario. Note that this number is NOT the same for all optimal scenarios.')
    parsesol_parser.set_defaults(func=dingII_parsesol.main)


def add_unimog_parsing_groups(parser):
    parser.add_argument('unimog', type=FileType('r'), help='The genomes provided in UniMOG format.')
    pairs = parser.add_mutually_exclusive_group()
    pairs.add_argument('-p', '--pair', type=str, nargs=2, help='Give the two names of the genomes you want to compare, as specified in the unimog header.')
    pairs.add_argument('-pn', '--pairnumber', type=int, nargs=2, help='Chose the two genomes via their position in the file (starting at 0). Default: 0,1')


def main():
    configure_logging()

    parser = ArgumentParser(description='ding II - an algorithm solving the genomic distance problem for '
                                                 'natural genomes, in which any marker may occur an arbitrary number'
                                                 ' of times')

    # Create subparsers
    subparsers = parser.add_subparsers(dest='subcommand', help='Available subcommands')

    # Define the generate subcommand
    define_generate_subcommand(subparsers)

    # Define the parsesol subcommand
    define_parsesol_subcommand(subparsers)

    args = parser.parse_args()

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    sys.exit(main())
