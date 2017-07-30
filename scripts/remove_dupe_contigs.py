import sys

if __name__ == "__main__":
    filename = sys.argv[1]
    with open(filename, 'r') as fp:
        records = {}
        for line in fp:
            key = line
            value = fp.next()
            records[key] = value.rstrip()
