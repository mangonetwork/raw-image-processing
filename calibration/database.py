import psycopg2
import sys
import os

# Connect to an existing database
conn = psycopg2.connect("dbname=test user=postgres")

# Open a cursor to perform database operations
cur = conn.cursor()


cur.close()
conn.close()