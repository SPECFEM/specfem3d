#!python
# Opening fault cracks
# Input : surface up and down.
import cubit 

class fault_input:
   def __init__(self,id,surface_u,surface_d): 
       self.id = id
       self.surface_u = surface_u
       self.surface_d = surface_d
       self.name = 'MESH/fault_file_'+str(id)+'.dat'

       quads_Aup,quads_Adp = save_cracks(self.name,self.surface_u,self.surface_d)
       #Unpacking list.
       quads_Au=unpack_list(quads_Aup)
       quads_Ad=unpack_list(quads_Adp)

       print 'len(Au):',len(quads_Au)
       print 'len(Ad):',len(quads_Ad)

       if not (len(quads_Au)==len(quads_Ad)):
           print 'Number of elements for each fauld side up and down do not concide'
           sys.exit('goodbye')
           
       save_elements_nodes(self.name,quads_Au,quads_Ad)


def save_cracks(name,list_surface_up,list_surface_down):
   quads_fault_up = [] 
   quads_fault_down = []
   for surface in list_surface_up   :
       quads_fault = cubit.get_surface_quads(surface)
       quads_fault_up.append(quads_fault)
   for surface in list_surface_down :
       quads_fault = cubit.get_surface_quads(surface)
       quads_fault_down.append(quads_fault)
    # TO DO : stop python properly in case fault nodes at both sides
    #         do not match.
    #   if len(quads_fault_u) != len(quads_fault_d): stop
    #
    # SAVING FAULT ELEMENTS AND NODES
   return quads_fault_up,quads_fault_down

def unpack_list(fault_list):
    list_fault = []
    for i in range(0,len(fault_list)):
        el=list(fault_list[i])
        for j in el:
            list_fault.append(j)
    return list_fault 
          
def save_elements_nodes(name,quads_fault_u,quads_fault_d):
   fault_file = open(name,'w')
   txt =''
   list_hex=cubit.parse_cubit_list('hex','all')
   txt='%10i %10i\n' % (len(quads_fault_u),len(quads_fault_d))
   fault_file.write(txt)
   
   dic_quads_fault_u = dict(zip(quads_fault_u,quads_fault_u)) 
   dic_quads_fault_d = dict(zip(quads_fault_d,quads_fault_d)) 
   
   # FAULT SIDE DOWN
   for h in list_hex: 
       faces = cubit.get_sub_elements('hex',h,2)  
       for f in faces:
           if dic_quads_fault_d.has_key(f): 
              nodes=cubit.get_connectivity('Face',f)
   #           print 'h,fault nodes side down :',h,nodes[0],nodes[1],nodes[2],nodes[3]
              txt='%10i %10i %10i %10i %10i\n' % (h,nodes[0],\
                                                nodes[1],nodes[2],nodes[3])
              fault_file.write(txt)

   # FAULT SIDE UP
   for h in list_hex: 
       faces = cubit.get_sub_elements('hex',h,2)  
       for f in faces:
           if dic_quads_fault_u.has_key(f): 
              nodes=cubit.get_connectivity('Face',f)
   #           print 'h,fault nodes side up :',h,nodes[0],nodes[1],nodes[2],nodes[3]
              txt='%10i %10i %10i %10i %10i\n' % (h,nodes[0],\
                                                nodes[1],nodes[2],nodes[3])
              fault_file.write(txt)

   fault_file.close()
