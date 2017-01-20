cubit.cmd('group \'nnn\' add Node in face in group lateral')
tol=5
empty=False
def mnode(tol,clean=True):
  print tol
  cubit.cmd("topology check coincident node node in group nnn tolerance "+str(tol)+" highlight brief result group 'chn'")
  group_exist=cubit.get_id_from_name("chn")
  if not group_exist:
    print 'no nodes in this tolerance range'
  else:
    merging_nodes=cubit.get_group_nodes(group_exist)
    cubit.cmd('draw group lateral')
    cubit.cmd('high group chn')
    print 'merging ',len(merging_nodes),' nodes.....'
    cubit.cmd("equivalence node in chn tolerance "+str(tol*2))
    if clean:
      cubit.cmd("group nnn remove node in group chn")
      cubit.cmd("delete Group chn")
    ilateral_nodes=cubit.get_id_from_name('nnn')
    lateral_nodes=cubit.get_group_nodes(ilateral_nodes)
    print len(lateral_nodes)
    if len(lateral_nodes) == 0:
      empty=True
  cubit.cmd('draw group lateral')
  cubit.cmd('high group nnn')
  cubit.cmd('quality hex all jacobian global high 0 draw mesh draw add')

for tol in range(1,16):
  mnode(tol)
