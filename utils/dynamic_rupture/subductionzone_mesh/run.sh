#!/bin/bash
python process_slab_rotate.py
trelis -nographics -input ./create_chunks_mesh.jou
python exportmesh.py
