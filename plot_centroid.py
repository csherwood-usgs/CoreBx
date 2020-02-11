# plot_centroids.py
# called by Process_CoreBx_multi_v3

cols=['#fdbe85','#fd8d3c','#e6550d','#a63603']
cbar = np.nanmean(c,1)
cbarall = np.nanmean(call,1)
fig,ax = plt.subplots(1)
plt.plot(cbar[:,0],cbar[:,1],'--',ms=12,c='dimgray')
#plt.plot(cbar[0,0],cbar[0,1],'ok',ms=12)
plt.plot(cbar[0,0],cbar[0,1],'o',ms=12,c=cols[0],label=dates[0])
#plt.plot(cbar[1,0],cbar[1,1],'ok',ms=12)
plt.plot(cbar[1,0],cbar[1,1],'o',ms=12,c=cols[1],label=dates[1])
#plt.plot(cbar[2,0],cbar[2,1],'ok',ms=12)
plt.plot(cbar[2,0],cbar[2,1],'o',ms=12,c=cols[2],label=dates[2])
#plt.plot(cbar[3,0],cbar[3,1],'ok',ms=12)
plt.plot(cbar[3,0],cbar[3,1],'o',ms=12,c=cols[3],label=dates[3])
#plt.gca().invert_xaxis()
plt.xlim((300.,150.))
plt.ylim((0., 2.5))
plt.ylabel('Elevation (m NAVD88)')
plt.xlabel('Cross-shore Distance (m)')
plt.legend()
plt.title('{} Centroids - {}'.format(r['name'],holes[ihole]))
fig_name = '{}_centroid_{}.svg'.format(r['name'],holes[ihole])
plt.savefig(fig_name,bbox_inches='tight', format='svg')
