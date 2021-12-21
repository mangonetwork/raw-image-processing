## Template obtained from Steven 02-01-2013
## Modified by Maria Rangel 02-01-203

from app import db
from app.users import constants as USER

class User(db.Model):

  __tablename__ = 'users_user'
  id = db.Column(db.Integer, primary_key = True)
  firstname = db.Column(db.String(20), unique = False)
  lastname = db.Column(db.String(30), unique = False)
  username = db.Column(db.String(20), unique = True)
  affiliation = db.Column(db.String(30), unique = False)
  email = db.Column(db.String(120), unique = True)
  password = db.Column(db.String(20), unique = False)
  role = db.Column(db.SmallInteger, default=USER.ADMIN)
  status = db.Column(db.SmallInteger, default=USER.NEW)

  def __init__(self, firstname=None, lastname=None, affiliation=None, username=None, email=None, password=None, role=None):
    self.firstname = firstname
    self.lastname = lastname
    self.username = username
    self.affiliation = affiliation
    self.email = email
    self.password = password
    self.role = role

  def getStatus(self):
    return USER.STATUS[self.status]

  def getRole(self):
    return USER.ROLE[self.role]

  def __repr__(self):
      return '<User %r>' % (self.username)
  
