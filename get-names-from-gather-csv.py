#! /usr/bin/env python
import csv
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('csv', nargs='+')
    parser.add_argument('-o', '--output')
    args = parser.parse_args()

    assert args.output, 'please specify an output with -o'

    all_names = set()
    for filename in args.csv:
        with open(filename, 'rt') as fp:
            rows = list(csv.DictReader(fp))

        if not rows or 'name' not in rows[0]:
            print('skipping', filename)
            continue

        names = [ r['name'] for r in rows ]
        all_names.update(names)

    print('writing {} names to "{}"'.format(len(all_names), args.output))

    with open(args.output, 'wt') as fp:
        fp.write("\n".join(all_names))


if __name__ == '__main__':
    main()
