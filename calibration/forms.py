#from flask.ext.wtf import Form, TextField, PasswordField, BooleanField, RecaptchaField
#from flask.ext.wtf import Required, Email, EqualTo
#from flask.ext.wtf import Form
from wtforms import TextField,PasswordField,BooleanField
from wtforms.validators import Required,Email,EqualTo

class LoginForm(Form):
    email = TextField('Email address', [Required(), Email()])
    password = PasswordField('Password', [Required()])

class RegisterForm(Form):
    firstname = TextField('First Name', [Required()])
    lastname = TextField('Last Name', [Required()])
    username = TextField('Username', [Required()])
    email = TextField('Email address', [Required(), Email()])
    password = PasswordField('Password', [Required()])
    confirm = PasswordField('Repeat Password', [Required(),EqualTo('password', message='Passwords must match')])
    affiliation = TextField('Affiliation', [Required()])
    accept_tos = BooleanField('I accept the TOS', [Required()])
    #recaptcha = RecaptchaField()

class LogoutForm(Form):
    pass

class ChangePassword(Form):
    current_password =  PasswordField('Current Password', [Required()])
    new_password = PasswordField('New Password', [Required()])
    confirm = PasswordField('Repeat Password', [Required(),EqualTo('new_password', message='Passwords must match')])
    
class AdminLoginForm(Form):
    email = TextField('Email address', [Required(), Email()])
    password = PasswordField('Password', [Required()])

class AdminRegisterForm(Form):
    firstname = TextField('First Name', [Required()])
    lastname = TextField('Last Name', [Required()])
    email = TextField('Email address', [Required(), Email()])
    password = PasswordField('Password', [Required()])
    confirm = PasswordField('Repeat Password', [Required(),EqualTo('confirm', message='Passwords must match')])
    accept_tos = BooleanField('I accept the TOS', [Required()])
    #recaptcha = RecaptchaField()

class AdminLogoutForm(Form):
    pass
