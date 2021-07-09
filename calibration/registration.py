from flask.ext.wtf import Form, TextField, PasswordField, BooleanField, RecaptchaField
from flask.ext.wtf import Required, Email, EqualTo
from flask.ext.login import LoginManager
from flask.ext.openid import OpenID
import os

class User(db.Model):
    id = db.Column(db.Integer, primary_key = True)
    nickname = db.Column(db.String(64), unique = True)
    email = db.Column(db.String(120), unique = True)
    role = db.Column(db.SmallInteger, default = ROLE_USER)
    posts = db.relationship('Post', backref = 'author', lazy = 'dynamic')

    def is_authenticated(self):
        return True

    def is_active(self):
        return True

    def is_anonymous(self):
        return False

    def get_id(self):
        return unicode(self.id)

    def __repr__(self):
        return '<User %r>' % (self.nickname)

class LoginForm(Form):
    email = TextField('Email address', [Required(), Email()])
    password = PasswordField('Password', [Required()])

class RegisterForm(Form):
    name = TextField('Name', [Required()])
    email = TextField('Email address', [Required(), Email()])
    password = PasswordField('Password', [Required()])
    confirm = PasswordField('Repeat Password', [Required(),EqualTo('confirm', message='Passwords must match')])
    accept_tos = BooleanField('I accept the TOS', [Required()])
    recaptcha = RecaptchaField()

class LogoutForm(Form):
    pass