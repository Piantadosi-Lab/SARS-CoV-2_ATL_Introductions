#!/usr/bin/env python3
# short wrapper script around numeric_from_datetime
# function which reads in datetime objects
# from stdin and prints numeric date to stdout
# takes a single argument to define date format
#echo "2020-03-31" | python3 scripts/numeric_from_datetime.py

def run():
    from treetime.utils import numeric_date
    import sys
    import argparse
    import datetime
    parser = argparse.ArgumentParser()
    parser.add_argument('--fmt', 
        default='%Y-%m-%d')
    args = parser.parse_args()
    for line in sys.stdin:
        line = line.strip()
        line_dt = datetime.datetime.strptime(line, args.fmt)
        converted_line = numeric_date(line_dt)
        print(converted_line, file=sys.stdout)




if __name__ == "__main__":
    run()
