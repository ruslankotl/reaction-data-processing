from datetime import datetime
from sys import argv

def main(files, destination):
    now = datetime.now()
    datm = now.strftime('$DATM %Y-%m-%d %H:%M:%S\n')
    header = '$RDFILE 1\n'+datm
    with open(destination, 'w') as dest_file:
        dest_file.write(header)
        for file in files:
            with open(file, 'r') as rdf:
                lines = iter(rdf.readlines())
                current = str()
                markush = False
                for line in lines:
                    if line.startswith('$RFMT'):
                        if not markush:
                            dest_file.write(current)
                        markush = False
                        current = str()
                    if 'MARKUSH' in line:
                        markush = True
                    current+=line
                dest_file.write(current)
    
    
if __name__ == "__main__":
    main(argv[1:-1], argv[-1])