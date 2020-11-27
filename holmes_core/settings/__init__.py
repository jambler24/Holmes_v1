from .base import *
import os
# you need to set "myproject = 'prod'" as an environment variable
# in your OS (on which your website is hosted)

'''
if os.environ['myproject'] == 'prod':
	from .prod import *
else:
	from .dev import *

'''

if os.environ.get('DJANGO_CONFIGURATION') == 'prod':
	from .prod import *
else:
	from .dev import *
