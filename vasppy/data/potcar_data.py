import yaml
import os

my_path = os.path.dirname(__file__)

potcar_md5sum_data = {}
potcar_nelect = {}

potcar_sets = [ 'PBE', 'PBE_52', 'PBE_54', 'GGA', 'USPP_GGA', 'LDA', 'LDA_52', 'LDA_54', 'USPP_LDA' ]

for potcar_set in potcar_sets:
    with open( '{}_md5.yaml'.format( os.path.join( my_path, potcar_set ) ), 'r' ) as stream:
        potcar_md5sum_data[ potcar_set ] = yaml.load( stream, Loader=yaml.SafeLoader  )

for potcar_set in potcar_sets:
    with open( '{}_nelect.yaml'.format( os.path.join( my_path, potcar_set ) ), 'r' ) as stream:
        potcar_nelect[ potcar_set ] = yaml.load( stream, Loader=yaml.SafeLoader  )
