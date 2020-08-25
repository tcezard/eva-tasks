# Context

https://www.ebi.ac.uk/panda/jira/browse/EVA-2124

We need to change chromosome names to chromosome accessions (and contigs).

The idea is similar to EVA-2068:
- Scan all documents in an assembly
- Filter by affected studies
- Drop the document and insert a similar document with the contig replaced to genbank, using assembly reports


# Set up

if you use a virtual environment:
```
virtualenv  -p python3.8 venv
source venv/bin/activate
pip install -r requirements.txt
```
you might need to install the next packages (specially for psycopg):
```
sudo apt-get install python3.8-dev postgresql-server-dev-10 postgresql-common
```


# Run

Without this, python won't find all the source files:
```
export PYTHONPATH=.
```
