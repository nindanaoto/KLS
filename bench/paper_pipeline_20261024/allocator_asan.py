from common import *
for role,version in [('recovery','v6'),('subtree','compact-v11')]:
 runner=Runner(role+'-'+version+'-asan-allocator','16-23');B=build(role,version,'asan')
 runner.run('build',['cmake','--build',B,'-j8','--target','kls_klu_allocator_ownership'])
 runner.run('test',[B/'kls_klu_allocator_ownership'])
 (O/(runner.name+'-completed.json')).write_text(json.dumps(dict(passed=True,scope='Four concurrent owners, tagged allocator families, allocation-failure rollback; ASan/UBSan/leaks',processes=runner.rows),indent=2))
