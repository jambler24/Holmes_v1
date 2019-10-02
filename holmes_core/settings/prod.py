
# Database
# https://docs.djangoproject.com/en/2.2/ref/settings/#databases

DATABASES = {
    'default': {
         'NAME': '/db/db.sqlite3',
         'ENGINE': 'django.db.backends.sqlite3',
     },

}


# Container Paths

ANNO_FOLDER = '/annotations/'
BAM_FILES_DIR = '/bam_files/'
REF_GENOME_FOLDER = '/ref_genome/'

