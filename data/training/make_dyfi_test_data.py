from mapping_tools import reckon


num_points = 4
point_dist = 50. # km
mmi_val = 9.5

angle_sep = 360. / num_points

csvtxt = 'LAT,LON,MMI\n'
for i in range(0, num_points):
    brngd = 0 + i*angle_sep
    ptlo, ptla = reckon(-26., 129., point_dist, brngd)
    
    csvtxt += ','.join((str('%0.4f' % ptla), str('%0.4f' % ptlo), str('%0.1f' % mmi_val))) + '\n'
    
f = open('test_data.csv', 'w')
f.write(csvtxt)
f.close()