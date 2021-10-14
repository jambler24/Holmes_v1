
# Database
# https://docs.djangoproject.com/en/2.2/ref/settings/#databases

DATABASES = {
    'default': {
         'NAME': 'db.sqlite3',
         'ENGINE': 'django.db.backends.sqlite3',
     },
}


'''
    'postgres': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': 'postgres',
        'USER': 'postgres',
        'HOST': 'db',
        'PORT': '5432',
    }
'''

# local dev paths

print('Dev mode on')

ANNO_FOLDER = '/Users/panix/Library/Mobile Documents/com~apple~CloudDocs/programs/Holmes/Docker_things/annotations/'
BAM_FILES_DIR = '/Users/panix/Library/Mobile Documents/com~apple~CloudDocs/programs/Holmes/Docker_things/bam_files/'
VARIANT_FOLDER = '/Users/panix/Library/Mobile Documents/com~apple~CloudDocs/programs/Holmes/Docker_things/variant_files/'

REF_GENOME_FOLDER = '/Users/panix/Library/Mobile Documents/com~apple~CloudDocs/programs/Holmes/Docker_things/ref_genome/'
#REF_GENOME_FOLDER = '/Volumes/External/genomes/hg19/'
