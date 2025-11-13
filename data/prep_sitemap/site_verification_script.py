#!/usr/bin/env python3
"""
Parse sitemap XYZ files and determine which matches docking SITE2
"""
import statistics

def parse_xyz_file(filename):
    """Parse XYZ file and return list of coordinates"""
    coords = []
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Skip first two lines (number of atoms and title)
    for line in lines[2:]:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) >= 4:
            try:
                # Format: ELEMENT X Y Z
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                coords.append((x, y, z))
            except:
                pass
    
    return coords

def calculate_centroid(coords):
    """Calculate centroid of a list of coordinates"""
    if not coords:
        return None
    
    x_avg = statistics.mean([c[0] for c in coords])
    y_avg = statistics.mean([c[1] for c in coords])
    z_avg = statistics.mean([c[2] for c in coords])
    
    return (x_avg, y_avg, z_avg)

def distance(coord1, coord2):
    """Calculate Euclidean distance between two 3D coordinates"""
    return ((coord1[0]-coord2[0])**2 + 
            (coord1[1]-coord2[1])**2 + 
            (coord1[2]-coord2[2])**2)**0.5

# Parse both site files
site2_coords = parse_xyz_file('/mnt/user-data/uploads/site2.xyz')
site3_coords = parse_xyz_file('/mnt/user-data/uploads/site3.xyz')

# Calculate centroids
site2_center = calculate_centroid(site2_coords)
site3_center = calculate_centroid(site3_coords)

# Docking coordinates
SITE1_DOCK = (27.991873, 18.851075, 15.488124)
SITE2_DOCK = (23.255162, 27.655522, 14.632348)

print("="*80)
print("SITEMAP SITE CENTER ANALYSIS")
print("="*80)
print()

print(f"sitemap_1_site_2:")
print(f"  Number of points: {len(site2_coords)}")
print(f"  Center (centroid): ({site2_center[0]:.3f}, {site2_center[1]:.3f}, {site2_center[2]:.3f})")
print()

print(f"sitemap_1_site_3:")
print(f"  Number of points: {len(site3_coords)}")
print(f"  Center (centroid): ({site3_center[0]:.3f}, {site3_center[1]:.3f}, {site3_center[2]:.3f})")
print()

print("="*80)
print("COMPARISON WITH DOCKING COORDINATES")
print("="*80)
print()

print(f"Docking SITE1: ({SITE1_DOCK[0]:.3f}, {SITE1_DOCK[1]:.3f}, {SITE1_DOCK[2]:.3f})")
print(f"Docking SITE2: ({SITE2_DOCK[0]:.3f}, {SITE2_DOCK[1]:.3f}, {SITE2_DOCK[2]:.3f})")
print()

# Calculate distances to SITE2
dist_site2_to_dock2 = distance(site2_center, SITE2_DOCK)
dist_site3_to_dock2 = distance(site3_center, SITE2_DOCK)

# Calculate distances to SITE1 (for reference)
dist_site2_to_dock1 = distance(site2_center, SITE1_DOCK)
dist_site3_to_dock1 = distance(site3_center, SITE1_DOCK)

print("Distance calculations:")
print("-"*80)
print(f"\nsitemap_1_site_2 to docking SITE1: {dist_site2_to_dock1:.3f} Å")
print(f"sitemap_1_site_2 to docking SITE2: {dist_site2_to_dock2:.3f} Å")

print(f"\nsitemap_1_site_3 to docking SITE1: {dist_site3_to_dock1:.3f} Å")
print(f"sitemap_1_site_3 to docking SITE2: {dist_site3_to_dock2:.3f} Å")

print()
print("="*80)
print("CONCLUSION")
print("="*80)
print()

if dist_site2_to_dock2 < dist_site3_to_dock2:
    print(f"✓ sitemap_1_site_2 MATCHES docking SITE2!")
    print(f"  Distance: {dist_site2_to_dock2:.3f} Å")
    print()
    print("  SiteMap site_2 → Your docking SITE2")
else:
    print(f"✓ sitemap_1_site_3 MATCHES docking SITE2!")
    print(f"  Distance: {dist_site3_to_dock2:.3f} Å")
    print()
    print("  SiteMap site_3 → Your docking SITE2")

print()
print("FINAL MAPPING:")
print("-"*80)
print("✓ sitemap_1_site_1 → Docking SITE1 (cryptic pocket)")
if dist_site2_to_dock2 < dist_site3_to_dock2:
    print("✓ sitemap_1_site_2 → Docking SITE2 (alternative site)")
    print("  sitemap_1_site_3 → Not used in docking")
else:
    print("  sitemap_1_site_2 → Not used in docking")
    print("✓ sitemap_1_site_3 → Docking SITE2 (alternative site)")
