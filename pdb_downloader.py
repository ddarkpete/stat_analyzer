import urllib2
from os import path



url = 'http://www.rcsb.org/pdb/rest/search'

query = open('search_query.xml').read()


print "query:\n", query

print "querying PDB...\n"

req = urllib2.Request(url, data=query)

f = urllib2.urlopen(req)

result = f.read()
ids = []

if result:

    print  result
    ids = result.split('\n')
    del ids[-1]
    print len(ids)
    for i in ids:
        pdb_url = 'https://files.rcsb.org/download/' + i  + '.pdb'
        try:
            pdb_file = urllib2.urlopen(pdb_url)
            save_file = 'pdbs/' + i + '.pdb'
            if path.exists(save_file):
                continue
            else:
                with open(save_file , 'w') as output:
                    output.write(pdb_file.read())
        except urllib2.HTTPError, e:
            if e.code == 404:
                print "No PDB file for id: " + i
                pass




else:

    print "Failed to retrieve results" 



