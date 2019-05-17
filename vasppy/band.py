import warnings

def handle_occupancy( occupancy, negative_occupancies='warn' ):
    valid_negative_occupancies = [ 'warn', 'raise', 'ignore', 'zero' ]
    if negative_occupancies not in valid_negative_occupancies:
        raise ValueError( "valid options for negative_occupancies are {}".format( valid_negative_occupancies ) )
    if occupancy < 0:
       if negative_occupancies == 'warn':
           warnings.warn( "One or more occupancies in your PROCAR file are negative." )
       elif negative_occupancies == 'raise':
           raise ValueError( "One or more occupancies in your PROCAR file are negative." )
       elif negative_occupancies == 'ignore':
           pass
       elif negative_occupancies == 'zero':
           occupancy = 0.0
    return occupancy

class Band():

    def __init__( self, index, energy, occupancy, negative_occupancies='warn' ):
        self.index = index
        self.energy = energy
        self.occupancy = handle_occupancy( occupancy, negative_occupancies=negative_occupancies )

    def __eq__( self, other ):
        return ( ( self.index == other.index ) &
                 ( self.energy == other.energy ) &
                 ( self.occupancy == other.occupancy ) )
