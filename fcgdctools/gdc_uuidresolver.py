"""
This module prototypes a GDC uuid resolver.

It takes as input a TSV file whose first column contains 
GDC file uuids and second column contains Google Cloud Storage
URLs.
"""

import dbm


class UuidResolver:

    """A GDC uuid resolver

    Attributes:
        tsvFile (str): path of TSV file that will be used to populate the
            the persistent key-value store (dbm).  Using dbm due to the potentially
            large numbers of keys needed.

        unknownResponse (str): string to return if uuid is not recognized
    """
    def __init__(self, tsvFile, unknownResponse):
        self.db_filename = "uuid_to_url"
        self.unknownResponse = unknownResponse
        db = dbm.open(self.db_filename, 'n')

        with open(tsvFile) as f:
            i = 0
            for line in f:
                i += 1
                uuid_url_pair = line.rstrip().split('\t')
                try:
                    db[uuid_url_pair[0]] = uuid_url_pair[1]
                except IndexError:
                    print(i)
                    print(uuid_url_pair)
                except Exception as e:
                    #print(i)
                    #print(uuid_url_pair)
                    #print("Unexpected error:", e)
                    # believe this is a threading issue in ndbm.  
                    # I found an on-line recommendation that if I sleep for 1 second
                    # and try again, would work.  It appears to
                    db[uuid_url_pair[0]] = uuid_url_pair[1]
            print("number of rows processed = {0}".format(i))
        db.close()

    def getURL(self, uuid):
        db = dbm.open(self.db_filename, 'r')
        try:
            url = db[uuid].decode("utf-8")
        except KeyError:
            url = self.unknownResponse
        db.close()
        return url
