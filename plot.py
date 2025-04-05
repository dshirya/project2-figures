import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from preprocess import parse_formula

def plot_scatter(plotting_df):
    """
    Creates a scatter plot of compound properties.
    
    Parameters:
        plotting_df (pd.DataFrame): DataFrame with columns:
            - "Formula"
            - "Neutron cross section"
            - "Neutron mass absorption"
    
    The x-axis corresponds to "Neutron cross section" and the y-axis to "Neutron mass absorption".
    Circle markers are used for the scatter plot.
    """
    # Create a high-resolution figure (300 dpi) without grid
    plt.figure(figsize=(8, 6), dpi=300)
    
    # Separate data for "Gd10RuCd3" and all other compounds
    gd_mask = plotting_df["Formula"] == "Gd10RuCd3"
    other_mask = ~gd_mask
    
    # Plot the other compounds in the default color
    plt.scatter(plotting_df.loc[other_mask, "Neutron cross section"],
                plotting_df.loc[other_mask, "Neutron mass absorption"],
                marker="o", alpha=0.2, color="black", edgecolors="none")
    
    # Plot Gd10RuCd3 in red
    # if gd_mask.any():
    #     plt.scatter(plotting_df.loc[gd_mask, "Neutron cross section"],
    #                 plotting_df.loc[gd_mask, "Neutron mass absorption"],
    #                 marker="o", color="red", alpha=0.2, edgecolors="none")
    
    plt.xlabel("Neutron Cross Section, barn", fontsize=14)
    plt.ylabel("Neutron Mass Absorption, m²/kg", fontsize=14)
    plt.tick_params(axis='both', labelsize=14)
  
    plt.savefig("plot.png", dpi=500, bbox_inches='tight')
    plt.show()

def plot_filtered_scatter(plotting_df):
    """
    Creates a scatter plot of compound properties after filtering.
    
    Only compounds with:
      - "Neutron cross section" > 31000 and
      - "Neutron mass absorption" > 5
    are plotted.
    
    Parameters:
        plotting_df (pd.DataFrame): DataFrame with columns:
            - "Formula"
            - "Neutron cross section"
            - "Neutron mass absorption"
    
    The x-axis corresponds to "Neutron cross section" and the y-axis to "Neutron mass absorption".
    Circle markers are used for the scatter plot.
    "Gd10RuCd3" is plotted in red, while all other compounds are plotted in the default blue.
    Each point is labeled with its formula.
    """
    # Filter the DataFrame based on the conditions
    filtered_df = plotting_df[(plotting_df["Neutron cross section"] > 31000) &
                                (plotting_df["Neutron mass absorption"] > 5)]
    
    # Create a high-resolution figure (300 dpi)
    plt.figure(figsize=(8, 6), dpi=300)
    
    # Separate data for "Gd10RuCd3" and all other compounds in the filtered DataFrame
    gd_mask = filtered_df["Formula"] == "Gd10RuCd3"
    other_mask = ~gd_mask
    
    # Plot the other compounds in default blue
    plt.scatter(filtered_df.loc[other_mask, "Neutron cross section"],
                filtered_df.loc[other_mask, "Neutron mass absorption"],
                marker="o", alpha=0.2, color="#0348a1", edgecolors="none")
    
    # Plot "Gd10RuCd3" in red if it exists
    # if gd_mask.any():
    #     plt.scatter(filtered_df.loc[gd_mask, "Neutron cross section"],
    #                 filtered_df.loc[gd_mask, "Neutron mass absorption"],
    #                 marker="o", color="red", alpha=0.2, edgecolors="none")
    
    # Add labels to each marker with an offset to avoid overlapping the marker
    for idx, row in filtered_df.iterrows():
        plt.annotate(
            row["Formula"],
            (row["Neutron cross section"], row["Neutron mass absorption"]),
            textcoords="offset points",
            xytext=(5, 5),
            ha="left",
            fontsize=8
        )
    
    plt.xlabel("Neutron Cross Section, barn")
    plt.ylabel("Neutron Mass Absorption, m²/kg")
    
    plt.savefig("plot_filtered.png", dpi=500, bbox_inches='tight')
    plt.show()
    
def plot_zoom_scatter(plotting_df):
    """
    Creates a scatter plot of compound properties after filtering.
    
    Only compounds with:
      - "Neutron cross section" > 31000 and
      - "Neutron mass absorption" > 5
    are plotted.
    
    Color conditions:
      - If the sum of the element counts (indices) in the formula equals 1 or 10,
        the marker is colored red.
      - If the sum equals 14, the marker is colored green.
      - Otherwise, the marker is colored blue (#0348a1).
    
    Markers not colored red (blue or green) are labeled with their formula.
    
    Parameters:
        plotting_df (pd.DataFrame): DataFrame with columns:
            - "Formula"
            - "Neutron cross section"
            - "Neutron mass absorption"
    
    The x-axis corresponds to "Neutron cross section, barn" and the y-axis to "Neutron mass absorption, m²/kg".
    Circle markers are used with alpha=0.2 and no edge colors.
    """
    # Filter the DataFrame based on the conditions
    filtered_df = plotting_df[(plotting_df["Neutron cross section"] > 34500) &
                                (plotting_df["Neutron mass absorption"] > 5)]
    
    # Create a high-resolution figure (300 dpi)
    plt.figure(figsize=(8, 6), dpi=300)
    
    # Iterate over each row in the filtered DataFrame to plot each marker with the color condition
    for idx, row in filtered_df.iterrows():
        formula = row["Formula"]
        # Parse the formula to get element counts
        element_counts = parse_formula(formula)
        sum_indices = sum(element_counts.values())
        
        # Set color based on conditions:
        # red if sum == 1 or 10, green if sum == 14, otherwise blue.
        if sum_indices in [1]:
            color = "red"
        elif sum_indices == 14:
            color = "green"
        else:
            color = "#0348a1"
        
        plt.scatter(row["Neutron cross section"],
                    row["Neutron mass absorption"],
                    marker="o", alpha=0.2, color=color, edgecolors="none")
        
        # Add label for compounds that are not red

        plt.annotate(formula,
                         (row["Neutron cross section"], row["Neutron mass absorption"]),
                         textcoords="offset points",
                         xytext=(5, 5),
                         ha="left", fontsize=8)
    
    plt.xlabel("Neutron Cross Section, barn")
    plt.ylabel("Neutron Mass Absorption, m²/kg")
    
    # Save the figure with high resolution
    plt.savefig("plot_zoom_nolabel.png", dpi=500, bbox_inches='tight')
    plt.show()

def plot_zoom_color_scatter(plotting_df):
    """
    Creates a scatter plot of compound properties after filtering.
    
    Only compounds with:
      - "Neutron cross section" > 34500 and
      - "Neutron mass absorption" > 5
    are plotted.
    
    The compounds are colored based on group membership:
    
      Group 1: Solid solution
          ["Gd", "Gd0.99Zr0.01", "Gd0.985Zr0.015", "Gd0.98In0.02", "Cd0.03Gd0.97",
           "Gd0.95Zr0.05", "Gd0.93Pr0.07", "Gd0.9Th0.1", "Cd0.15Gd0.85",
           "Cd0.17Gd0.83", "Gd0.8Nd0.2", "Gd0.8Pr0.2", "Ce0.2Gd0.8",
           "Gd0.78La0.22", "Ce0.28Gd0.72"]
          Color: "#c3121e"  (Sangre)
      
      Group 2: Element inclusion
          ["C0.03Gd", "C0.33Gd", "C0.4Gd"]
          Color: "#0348a1"  (Neptune)
      
      Group 3: GdM3 binary
          ["Gd3Ir", "Gd3In", "Fe0.04Gd3In0.96", "Gd3In0.95Mn0.05",
           "Gd3In0.95Ni0.05", "Al0.08Gd3In0.92", "Gd3Rh", "CoGd3",
           "Co0.844Gd3V0.156", "Co0.5Gd3Ni0.5", "Co0.5Gd3Ru0.5",
           "Gd3Os", "Gd3Pt", "Co0.2Gd3Ru0.8", "Gd3Pd", "Gd3Ni",
           "Gd3Ru", "AlGd3"]
          Color: "#ffb01c"  (Pumpkin)
      
      Group 4: Ternary
          ["Co2Gd14In3", "Cd3Gd10Os", "Gd10RuCd3", "Co3Gd14In2.7", "Fe1.29Gd5In0.79"]
          Color: "#027608"  (Clover)
      
      Group 5: Other binary
          ["CoGd9", "AlC0.1Gd3", "Gd5Ir2", "Gd5Ru2"]
          Color: "#1dace6"  (Cerulean)
    
    Compounds not found in any group are plotted in gray.
    
    Parameters:
        plotting_df (pd.DataFrame): DataFrame with columns:
            - "Formula"
            - "Neutron cross section"
            - "Neutron mass absorption"
    
    The x-axis corresponds to "Neutron cross section, barn" and the y-axis to "Neutron mass absorption, m²/kg".
    Circle markers are used with alpha=0.2 and no edge colors.
    """
    # Define groups as sets for faster membership testing
    group1 = {"Gd", "Gd0.99Zr0.01", "Gd0.985Zr0.015", "Gd0.98In0.02", "Cd0.03Gd0.97",
              "Gd0.95Zr0.05", "Gd0.93Pr0.07", "Gd0.9Th0.1", "Cd0.15Gd0.85",
              "Cd0.17Gd0.83", "Gd0.8Nd0.2", "Gd0.8Pr0.2", "Ce0.2Gd0.8",
              "Gd0.78La22", "Ce0.28Gd0.72"}  # Note: Check spelling consistency (e.g., "Gd0.78La22" vs "Gd0.78La0.22")
    # Adjust group1 if needed. For now, I'll assume the provided list is correct.
    group1 = {"Gd", "Gd0.99Zr0.01", "Gd0.985Zr0.015", "Gd0.98In0.02", "Cd0.03Gd0.97",
              "Gd0.95Zr0.05", "Gd0.93Pr0.07", "Gd0.9Th0.1", "Cd0.15Gd0.85",
              "Cd0.17Gd0.83", "Gd0.8Nd0.2", "Gd0.8Pr0.2", "Ce0.2Gd0.8",
              "Gd0.78La0.22", "Ce0.28Gd0.72"}
    
    group2 = {"C0.03Gd", "C0.33Gd", "C0.4Gd"}
    
    group3 = {"Gd3Ir", "Gd3In", "Fe0.04Gd3In0.96", "Gd3In0.95Mn0.05",
              "Gd3In0.95Ni0.05", "Al0.08Gd3In0.92", "Gd3Rh", "CoGd3",
              "Co0.844Gd3V0.156", "Co0.5Gd3Ni0.5", "Co0.5Gd3Ru0.5",
              "Gd3Os", "Gd3Pt", "Co0.2Gd3Ru0.8", "Gd3Pd", "Gd3Ni",
              "Gd3Ru", "AlGd3"}
    
    group4 = {"Co2Gd14In3", "Cd3Gd10Os", "Gd10RuCd3", "Co3Gd14In2.7", "Fe1.29Gd5In0.79"}
    
    group5 = {"CoGd9", "AlC0.1Gd3", "Gd5Ir2", "Gd5Ru2"}
    
    # Color mapping for groups:
    color_map = {
        1: "#c3121e",  # Sangre
        2: "#0348a1",  # Neptune
        3: "#ffb01c",  # Pumpkin
        4: "#027608",  # Clover
        5: "#1dace6"   # Cerulean
    }
    
    # Filter the DataFrame based on the conditions
    filtered_df = plotting_df[(plotting_df["Neutron cross section"] > 34500) & (plotting_df["Neutron cross section"] < 50500) &
                                (plotting_df["Neutron mass absorption"] > 5)]
    
    # Create a high-resolution figure (300 dpi)
    plt.figure(figsize=(8, 6), dpi=300)
    
    # Iterate over each row in the filtered DataFrame
    for idx, row in filtered_df.iterrows():
        formula = row["Formula"]
        # Determine group based on membership
        if formula in group1:
            color = color_map[1]
        elif formula in group2:
            color = color_map[2]
        elif formula in group3:
            color = color_map[3]
        elif formula in group4:
            color = color_map[4]
        elif formula in group5:
            color = color_map[5]
        else:
            color = "gray"  # default color for compounds not in any group
        
        plt.scatter(row["Neutron cross section"],
                    row["Neutron mass absorption"],
                    marker="o", alpha=0.5, color=color, edgecolors="none", s=130)

    # plt.xlabel("Neutron Cross Section, barn")
    # plt.ylabel("Neutron Mass Absorption, m²/kg")
    # legend_handles = [
    #     mpatches.Patch(color=color_map[1], label="Solid solution"),
    #     mpatches.Patch(color=color_map[2], label="Element inclusion"),
    #     mpatches.Patch(color=color_map[3], label="GdM3 binary"),
    #     mpatches.Patch(color=color_map[4], label="Ternary"),
    #     mpatches.Patch(color=color_map[5], label="Other binary")
    # ]
    # plt.legend(handles=legend_handles, loc="best")
    plt.tick_params(axis='both', labelsize=14)

    plt.savefig("plot_grouped.png", dpi=500, bbox_inches="tight")
    plt.show()