import bpy
import numpy as np
import os
import time
import json

#
# IMPORTANT
#
# Do not run this script natively. This script is meant to be run in Blender
# via one of the call routines
#

with open(os.path.join(os.path.dirname(__file__), 'manifest.json')) as f:
    manifest = json.load(f)

def main():
    # set the scene
    settings = {
        'resolution': 512,
        'camera_location': (-10,0,0),
        'camera_rotation': (np.pi/2,0,-np.pi/2),
        'camera_scale' : 10
    }
    set_environment(settings)

    # read molecule file and load it
    mol = read_xyz(manifest['xyzfile'])
    create_atoms(mol)
    create_bonds(mol)

    add_isosurface(manifest['mo_name'],
                   manifest['mo_neg_path'],
                   manifest['mo_pos_path'])

    render_scene(manifest['png_output'])

def add_isosurface(label, filename_neg, filename_pos):
    """
    Add the two isosurfaces
    """
    bpy.ops.import_mesh.ply(
        filepath=filename_neg
    )
    obj = bpy.context.object
    bpy.ops.object.shade_smooth()
    obj.data.materials.append(create_material('matneg', manifest['mo_colors']['neg'], alpha=0.5))
    obj.name = 'isosurface ' + label + '_neg'

    obj = bpy.ops.import_mesh.ply(
        filepath=filename_pos
    )
    bpy.ops.object.shade_smooth()
    obj = bpy.context.object
    obj.data.materials.append(create_material('matpos', manifest['mo_colors']['pos'], alpha=0.5))
    obj.name = 'isosurface ' + label + '_pos'

def create_atoms(mol):
    """
    Create atoms
    """
    for i,at in enumerate(mol):
        scale = manifest['atom_radii'][at[0]]
        bpy.ops.surface.primitive_nurbs_surface_sphere_add(
            radius=scale,
            enter_editmode=False,
            align='WORLD',
            location=at[1])
        obj = bpy.context.view_layer.objects.active
        obj.name = "atom-%s-%03i" % (at[0],i)
        bpy.ops.object.shade_smooth()

        # set a material
        mat = create_material(at[0], manifest['atom_colors'][at[0]])
        print(mat)
        obj.data.materials.append(mat)

def create_bonds(mol):
    """
    Create bonds between atoms
    """
    # set default orientation of bonds (fixed!)
    z = np.array([0,0,1])

    # add new bonds material if it does not yet exist
    matbond = create_material('bond', '222222')

    for i,at1 in enumerate(mol):
        r1 = np.array(at1[1])
        for j,at2 in enumerate(mol[i+1:]):
            r2 = np.array(at2[1])
            dist = np.linalg.norm(r2 - r1)

            # only create a bond if the distance is less than 1.5 A
            if dist < 2.5:
                axis = np.cross(z,r2-r1)
                if np.linalg.norm(axis) < 1e-5:
                    axis = np.array([0,0,1])
                    angle = 0.0
                else:
                    axis /= np.linalg.norm(axis)
                    angle = np.arccos(np.dot(r2-r1,z)/dist)

                bpy.ops.surface.primitive_nurbs_surface_cylinder_add(
                    enter_editmode=False,
                    align='WORLD',
                    location=tuple((r1 + r2) / 2)
                )

                obj = bpy.context.view_layer.objects.active
                obj.scale = (manifest['bond_thickness'], manifest['bond_thickness'], dist/2)
                obj.rotation_mode = 'AXIS_ANGLE'
                obj.rotation_axis_angle = (angle, axis[0], axis[1], axis[2])

                obj.name = "bond-%s-%03i-%s-%03i" % (at1[0],i,at2[0],j)
                bpy.ops.object.shade_smooth()
                obj.data.materials.append(matbond)

def set_environment(settings):
    """
    Specify canvas size, remove default objects, reset positions of
    camera and light, define film and set background
    """
    print('Set render engine to: CYCLES')
    bpy.context.scene.render.engine = 'CYCLES'
    print('Set rendering device to GPU')
    bpy.context.scene.cycles.device = 'GPU'
    bpy.context.scene.render.resolution_x = settings['resolution']
    bpy.context.scene.render.resolution_y = settings['resolution']
    print('Setting resolution to: ', settings['resolution'])
    bpy.context.scene.cycles.samples = 1024
    bpy.context.scene.cycles.tile_size = 2048

    # remove cube
    if 'Cube' in bpy.data.objects:
        o = bpy.data.objects['Cube']
        bpy.data.objects.remove(o, do_unlink=True)

    # set camera into default position
    bpy.data.objects['Camera'].location = tuple(settings['camera_location'])
    bpy.data.objects['Camera'].rotation_euler = tuple(settings['camera_rotation'])
    bpy.data.objects['Camera'].data.clip_end = 1000
    bpy.data.objects['Camera'].data.type = 'ORTHO'
    bpy.data.objects['Camera'].data.ortho_scale = settings['camera_scale']

    # set lights
    bpy.data.objects['Light'].data.type = 'AREA'
    bpy.data.objects['Light'].data.energy = 1e4
    bpy.data.objects['Light'].location = (-10,10,10)
    bpy.data.objects['Light'].rotation_euler = tuple(np.radians([55, 0, 225]))
    bpy.data.objects['Light'].data.shape = 'DISK'
    bpy.data.objects['Light'].data.size = 10

    # set film
    bpy.context.scene.render.film_transparent = True

    # set background
    bpy.data.worlds["World"].node_tree.nodes["Background"].inputs[0].default_value = (1,1,1,1)
    bpy.data.worlds["World"].node_tree.nodes["Background"].inputs[1].default_value = 1

def create_material(name, color, alpha=1.0):
    """
    Build a new material
    """
    # early exit if material already exists
    if name in bpy.data.materials:
        return bpy.data.materials[name]

    mat = bpy.data.materials.new(name)
    mat.use_nodes = True

    # set base color
    mat.node_tree.nodes["Principled BSDF"].inputs[0].default_value = hex2rgbtuple(color)

    # subsurface modifier
    mat.node_tree.nodes["Principled BSDF"].inputs[1].default_value = 0.2

    # set subsurface radii
    mat.node_tree.nodes["Principled BSDF"].inputs[2].default_value = (0.3,0.3,0.3)

    # set subsurface color
    mat.node_tree.nodes["Principled BSDF"].inputs[3].default_value = hex2rgbtuple('000000')

    # metallic
    mat.node_tree.nodes["Principled BSDF"].inputs[4].default_value = 0.3

    # roughness
    mat.node_tree.nodes["Principled BSDF"].inputs[7].default_value = 0.05

    # alpha
    mat.node_tree.nodes["Principled BSDF"].inputs[21].default_value = alpha

    return mat

def render_scene(outputfile, samples=512):
    """
    Render the scene
    """
    bpy.context.scene.cycles.samples = samples

    print('Start render')
    start = time.time()
    bpy.data.scenes['Scene'].render.filepath = outputfile
    bpy.ops.render.render(write_still=True)
    end = time.time()
    print('Finished rendering frame in %.1f seconds' % (end - start))

def read_xyz(filename):
    f = open(filename)
    nratoms = int(f.readline())
    f.readline() # skip line

    angtobohr = 1.8897259886

    mol = []
    for i in range(0,nratoms):
        pieces = f.readline().split()
        mol.append(
            [pieces[0],
            (
                float(pieces[1]) * angtobohr,
                float(pieces[2]) * angtobohr,
                float(pieces[3]) * angtobohr
            )]
        )

    return mol

def hex2rgbtuple(hexcode):
    """
    Convert 6-digit color hexcode to a tuple of floats
    """
    hexcode += "FF"
    hextuple = tuple([int(hexcode[i:i+2], 16)/255.0 for i in [0,2,4,6]])

    return tuple([color_srgb_to_scene_linear(c) for c in hextuple])

def color_srgb_to_scene_linear(c):
    """
    Convert RGB to sRGB
    """
    if c < 0.04045:
        return 0.0 if c < 0.0 else c * (1.0 / 12.92)
    else:
        return ((c + 0.055) * (1.0 / 1.055)) ** 2.4

if __name__ == '__main__':
    main()
