#!/usr/bin/env python3
# short wrapper script around numeric_from_datetime
# function which reads in datetime objects
# from stdin and prints numeric date to stdout
# takes a single argument to define date format


def run():
    from treetime.utils import datetime_from_numeric
    import sys
    import argparse
    import datetime
    #print(sys.stdin)
    parser = argparse.ArgumentParser()
    parser.add_argument('--fmt', 
        default='%Y-%m-%d')
    args = parser.parse_args()
    for line in sys.stdin:
        line = float(line.strip())
        line_dt = datetime_from_numeric(line)
        print(line_dt.strftime(args.fmt), file=sys.stdout)




if __name__ == "__main__":
    run()
