import agwsvalid

dgnf = agwsvalid.validator('dgnf')
m3   = agwsvalid.validator('m3')
dgwf = agwsvalid.validator('dgwf')

dgnfresult = dgnf.check(0 , .12, -.12, 0, 0, -.12, .12, 0)

m3result = m3.check(0 , .12, -.12, 0, 0, -.12, .12, 0)

dgwfresult = dgwf.check(0 , .12, -.12, 0, 0, -.12, .12, 0)


print(dgnfresult, m3result, dgwfresult)


