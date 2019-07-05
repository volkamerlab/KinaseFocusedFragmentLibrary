import apsw


# the extension is usually loaded right after the connection to the
# database
connection = apsw.Connection('chembldb.sql')
connection.enableloadextension(True)
connection.loadextension('~/chemicalite/')
connection.enableloadextension(False)
