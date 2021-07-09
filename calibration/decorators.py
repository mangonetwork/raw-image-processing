from functools import wraps

from flask import g, flash, redirect, url_for, request

from app.users import constants as USER

def requires_login(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):
        if g.user is None:
            flash(u'You need to be signed in for this page.')
            return redirect(url_for('users.login', next=request.path))
        return f(*args, **kwargs)
    return decorated_function

def adminrequires_login(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):
        if g.user is None:
            flash(u'You need to be signed in for this page.')
            return redirect(url_for('users.MANGO', next=request.path))
        if g.user.role != USER.ADMIN:
            flash(u'You do not have administrator priveledges to view this page.')
            return redirect(url_for('users.MANGO', next=request.path))
        return f(*args, **kwargs)
    return decorated_function
