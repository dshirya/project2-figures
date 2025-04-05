import pandas as pd
import matplotlib.pyplot as plt
import re

def read_compounds(file_path):
    df = pd.read_excel(file_path)
    return df

def parse_formula(formula):
    pattern = r'([A-Z][a-z]*)(\d*\.?\d*)'
    matches = re.findall(pattern, formula)
    elements = {}
    for elem, count in matches:
        count = float(count) if count else 1.0
        elements[elem] = elements.get(elem, 0) + count
    return elements

def read_elemental_properties(file_path):
    df = pd.read_excel(file_path)
    columns_to_keep = ["Symbol", "Atomic weight", "Neutron Cross Section", "Neutron Mass Absorption"]
    df = df[columns_to_keep]
    return df

def is_valid_formula(formula, valid_elements):
    parsed_elements = parse_formula(formula)
    return set(parsed_elements.keys()).issubset(valid_elements)

def filter_compounds(compounds_df, elemental_properties_df):
    # Create a set of valid element symbols from the elemental properties DataFrame
    valid_elements = set(elemental_properties_df["Symbol"].unique())
    
    # Apply the validity check for each formula in the first column
    mask = compounds_df.iloc[:, 0].apply(lambda formula: is_valid_formula(formula, valid_elements))
    filtered_df = compounds_df[mask].reset_index(drop=True)
    return filtered_df

def calculate_compound_neutron_cross_section(formula, elemental_properties_df):
    """
    Calculates the compound neutron cross section as a weighted average using 
    atomic percentages (based on the indices in the formula).
    
    Parameters:
        formula (str): The chemical formula (e.g., "H2O", "C6H12O6").
        elemental_properties_df (pd.DataFrame): DataFrame containing at least the 
            "Symbol" and "Neutron Cross Section" columns.
    
    Returns:
        tuple: (compound_cross_section, contributions)
            - compound_cross_section (float): The overall neutron cross section of the compound.
            - contributions (dict): A dictionary with each element's contribution, computed as
              atomic fraction multiplied by the element's neutron cross section.
    """
    parsed_elements = parse_formula(formula)
    total_atoms = sum(parsed_elements.values())
    # Calculate atomic fractions
    atomic_fractions = {element: count / total_atoms for element, count in parsed_elements.items()}
    
    # Retrieve neutron cross sections by element
    cross_sections = elemental_properties_df.set_index("Symbol")["Neutron Cross Section"].to_dict()
    
    contributions = {}
    compound_cross_section = 0.0
    for element, fraction in atomic_fractions.items():
        element_cross_section = cross_sections.get(element, 0)
        contributions[element] = fraction * element_cross_section
        compound_cross_section += contributions[element]
    
    return compound_cross_section, contributions

def calculate_compound_neutron_mass_absorption(formula, elemental_properties_df):
    """
    Calculates the compound neutron mass absorption using mass percentages.
    
    The mass percentages are computed based on the atomic weights of the elements (from the 
    elemental properties DataFrame). The compound neutron mass absorption is the weighted sum
    of each element's neutron mass absorption by its mass fraction.
    
    Parameters:
        formula (str): The chemical formula (e.g., "H2O", "C6H12O6").
        elemental_properties_df (pd.DataFrame): DataFrame containing at least the 
            "Symbol", "Atomic weight", and "Neutron Mass Absorption" columns.
    
    Returns:
        tuple: (compound_mass_absorption, contributions)
            - compound_mass_absorption (float): The overall neutron mass absorption of the compound.
            - contributions (dict): A dictionary with each element's contribution, computed as
              mass fraction multiplied by the element's neutron mass absorption.
    """
    parsed_elements = parse_formula(formula)
    
    # Get atomic weights and neutron mass absorptions
    atomic_weights = elemental_properties_df.set_index("Symbol")["Atomic weight"].to_dict()
    mass_absorptions = elemental_properties_df.set_index("Symbol")["Neutron Mass Absorption"].to_dict()
    
    # Calculate mass contributions of each element
    mass_contributions = {}
    total_mass = 0.0
    for element, count in parsed_elements.items():
        weight = atomic_weights.get(element, 0)
        mass = count * weight
        mass_contributions[element] = mass
        total_mass += mass
    
    # Compute mass fractions
    mass_fractions = {element: mass / total_mass for element, mass in mass_contributions.items()}
    
    contributions = {}
    compound_mass_absorption = 0.0
    for element, fraction in mass_fractions.items():
        element_absorption = mass_absorptions.get(element, 0)
        contributions[element] = fraction * element_absorption
        compound_mass_absorption += contributions[element]
    
    return compound_mass_absorption, contributions

def create_plotting_dataframe(compounds_df, elemental_properties_df):
    """
    Creates a DataFrame for plotting with columns: 
    'Formula', 'Neutron cross section', and 'Neutron mass absorption'.
    
    It calculates the compound neutron cross section and mass absorption for each valid formula.
    
    Parameters:
        compounds_df (pd.DataFrame): DataFrame containing chemical formulas in the first column.
        elemental_properties_df (pd.DataFrame): DataFrame with elemental properties.
    
    Returns:
        pd.DataFrame: A DataFrame with the required columns.
    """
    data = []
    # Set of valid elements from the elemental properties DataFrame.
    valid_elements = set(elemental_properties_df["Symbol"].unique())
    
    for formula in compounds_df.iloc[:, 0]:
        if is_valid_formula(formula, valid_elements):
            cross_section, _ = calculate_compound_neutron_cross_section(formula, elemental_properties_df)
            mass_absorption, _ = calculate_compound_neutron_mass_absorption(formula, elemental_properties_df)
            data.append({
                "Formula": formula,
                "Neutron cross section": cross_section,
                "Neutron mass absorption": mass_absorption
            })
    
    return pd.DataFrame(data)
