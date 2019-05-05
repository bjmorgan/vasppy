import yaml
import re
from collections import Counter

class Calculation:
    """
    class describing a single VASP calculation
    """
    
    def __init__( self, title, energy, stoichiometry ):
        """
        Initialise a Calculation object

        Args:
            title (Str): The title string for this calculation.
            energy (Float): Final energy in eV.
            stoichiometry (Dict{Str:Int}): A dict desribing the calculation stoichiometry,
                e.g. { 'Ti': 1, 'O': 2 }

        Returns:
            None
        """
        self.title = title
        self.energy = energy
        self.stoichiometry = Counter( stoichiometry )
        
    def __mul__( self, scaling ):
        """
        "Multiply" this Calculation by a scaling factor.
        Returns a new Calculation with the same title, but scaled energy and stoichiometry.

        Args:
            scaling (float): The scaling factor.

        Returns:
            (vasppy.Calculation): The scaled Calculation.
        """
        new_calculation = Calculation( title=self.title, energy=self.energy*scaling, stoichiometry=self.scale_stoichiometry( scaling ) )
        return new_calculation
        
    def __truediv__( self, scaling ):
        """
        Implements division by a scaling factor.
        Returns a new Calculation with the same title, but scaled energy and stoichiometry.

        Args:
            scaling (float): The scaling factor.

        Returns:
            (vasppy.Calculation): The scaled Calculation.
        """
        return self * ( 1 / scaling )
        
    def scale_stoichiometry( self, scaling ):
        """
        Scale the Calculation stoichiometry
        Returns the stoichiometry, scaled by the argument scaling.

        Args:
            scaling (float): The scaling factor.

        Returns:
            (Counter(Str:Int)): The scaled stoichiometry as a :obj:`Counter` of ``label: stoichiometry`` pairs
        """ 
        return { k:v*scaling for k,v in self.stoichiometry.items() }
    
def delta_E( reactants, products, check_balance=True ):
    """
    Calculate the change in energy for reactants --> products.
    
    Args:
        reactants (list(:obj:`vasppy.Calculation`)): A `list` of :obj:`vasppy.Calculation` objects. The initial state.
        products  (list(:obj:`vasppy.Calculation`)): A `list` of :obj:`vasppy.Calculation` objects. The final state.
        check_balance (:obj:`bool`, optional): Check that the reaction stoichiometry is balanced. Default is ``True``.

    Returns:
        (float) The change in energy.
    """
    if check_balance:
        if delta_stoichiometry( reactants, products ) != {}:
            raise ValueError( "reaction is not balanced: {}".format( delta_stoichiometry( reactants, products) ) )
    return sum( [ r.energy for r in products ] ) - sum( [ r.energy for r in reactants ] )

def delta_stoichiometry( reactants, products ):
    """
    Calculate the change in stoichiometry for reactants --> products.

    Args:
        reactants (list(:obj:`vasppy.Calculation`): A `list` of :obj:`vasppy.Calculation objects.` The initial state.
        products  (list(:obj:`vasppy.Calculation`): A `list` of :obj:`vasppy.Calculation objects.` The final state.

    Returns:
        (Counter): The change in stoichiometry.
    """ 
    totals = Counter()
    for r in reactants:
        totals.update( ( r * -1.0 ).stoichiometry )
    for p in products:
        totals.update( p.stoichiometry )
    to_return = {}
    for c in totals:
        if totals[c] != 0:
            to_return[c] = totals[c]
    return to_return

def energy_string_to_float( string ):
    """
    Convert a string of a calculation energy, e.g. '-1.2345 eV' to a float.

    Args:
        string (str): The string to convert.
  
    Return
        (float) 
    """
    energy_re = re.compile( "(-?\d+\.\d+)" )
    return float( energy_re.match( string ).group(0) )
    
def import_calculations_from_file( filename ):
    """
    Construct a list of :obj:`Calculation` objects by reading a YAML file.
    Each YAML document should include ``title``, ``stoichiometry``, and ``energy`` fields, e.g.::

        title: my calculation
        stoichiometry:
            - A: 1
            - B: 2
        energy: -0.1234 eV

    Separate calculations should be distinct YAML documents, separated by `---`
    
    Args:
        filename (str): Name of the YAML file to read.

    Returns:
        (dict(vasppy.Calculation)): A dictionary of :obj:`Calculation` objects. For each :obj:`Calculation` object, the ``title`` field from the YAML input is used as the dictionary key.
    """
    calcs = {}
    with open( filename, 'r' ) as stream:
        docs = yaml.load_all( stream, Loader=yaml.SafeLoader )
        for d in docs:
            stoichiometry = Counter()
            for s in d['stoichiometry']:
                stoichiometry.update( s )
            calcs[ d['title'] ] = Calculation( title=d['title'], 
                                               stoichiometry=stoichiometry, 
                                               energy=energy_string_to_float( d['energy'] ) )
    return calcs
