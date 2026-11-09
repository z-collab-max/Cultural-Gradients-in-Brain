"""Signed-median replacement curves and decile random-order tests."""
import numpy as np

def q(x):return dict(zip(('lower','median','upper'),map(float,np.quantile(x,[.025,.5,.975]))))
def recover(E,T):
    E=E-E.mean(1,keepdims=True);T=np.atleast_2d(T);T=T-T.mean(1,keepdims=True)
    ss=np.sum(T*T,1);cross=E@T.T
    return (ss[None,:]+cross)/np.sqrt(ss[None,:]*(np.sum(E*E,1)[:,None]+ss[None,:]+2*cross))


def replacement(samples,test_orders=2000):
    records=[];data={};rng=np.random.default_rng(2026100201)
    for k in (1,2):
        d=samples[k-1]; draws,m=d['beta_null'].shape
        theta=d['truth']/d['sd_um'];E=d['beta_null']/d['sd_um'];se=d['source_se']/d['sd_um']
        med=float(np.median(theta));top=np.argsort(-np.abs(theta/se),kind='stable')
        controls=np.array([rng.permutation(m) for _ in range(100)]);curves=[]
        for pct in range(100):
            n=round(m*pct/100)
            if n>=m:continue
            orders=np.vstack([top,controls]);T=np.tile(theta,(len(orders),1))
            for j,o in enumerate(orders):T[j,o[:n]]=med
            r=recover(E,T)
            curves.append({'percent':pct,'count':n,'ranked':q(r[:,0]),'random':q(r[:,1:]),'gap':q(r[:,1:].mean(1)-r[:,0])})
        # Random-order test at deciles.
        pct=np.arange(10,100,10);counts=np.rint(m*pct/100).astype(int);B=test_orders
        if np.any((counts<1)|(counts>=m)):raise ValueError("At least six distinct regional targets are required")
        orders=np.array([rng.permutation(m) for _ in range(B)]);delta=med-theta
        ec=E-E.mean(1,keepdims=True);ecross=ec@theta;enorm=np.sum(ec*ec,1);total=theta.sum();sq=np.sum(theta*theta)
        ref=np.empty((B,len(pct)))
        for start in range(0,B,20):
            o=orders[start:start+20];dd=delta[o]
            s1=total+np.cumsum(dd,1);s2=sq+np.cumsum(2*theta[o]*dd+dd*dd,1);S=(s2-s1*s1/m)[:,counts-1]
            cross=ecross[:,None,None]+np.cumsum(ec[:,o]*dd[None,:,:],2)[:,:,counts-1]
            rr=(S[None,:,:]+cross)/np.sqrt(S[None,:,:]*(S[None,:,:]+2*cross+enorm[:,None,None]))
            ref[start:start+len(o)]=rr.mean(0)
        T=np.tile(theta,(len(pct),1))
        for j,n in enumerate(counts):T[j,top[:n]]=med
        topmean=recover(E,T).mean(0);p=np.mean(ref<=topmean[None,:],axis=0)
        tests=[{'PC':k,'percent':int(x),'count':int(n),'top_mean_r':float(topmean[j]),'random_mean_r':float(ref[:,j].mean()),'one_sided_p':float(p[j])} for j,(x,n) in enumerate(zip(pct,counts))]
        data.update({f'pc{k}_null':ref,f'pc{k}_observed':topmean,f'pc{k}_test_pct':pct})
        records.append({'PC':k,'curves':curves,'tests':tests})

    # Holm across all 18 comparisons.
    flat=[z for r in records for z in r['tests']];raw=np.array([z['one_sided_p'] for z in flat]);order=np.argsort(raw);adj=np.empty(len(raw));running=0
    for rank,index in enumerate(order):running=max(running,min(1,(len(raw)-rank)*raw[index]));adj[index]=running
    for z,a in zip(flat,adj):z['holm_p']=float(a);z['holm_significant']=bool(a<.05)

    for k in (1,2):
        curves=records[k-1]['curves'];x=np.array([z['percent'] for z in curves])
        tl=np.array([z['ranked']['lower'] for z in curves]);tm=np.array([z['ranked']['median'] for z in curves]);tu=np.array([z['ranked']['upper'] for z in curves])
        rl=np.array([z['random']['lower'] for z in curves]);rm=np.array([z['random']['median'] for z in curves]);ru=np.array([z['random']['upper'] for z in curves])
        sig=np.array([z['percent'] for z in records[k-1]['tests'] if z['holm_significant']],float)
        ns=np.array([z['percent'] for z in records[k-1]['tests'] if not z['holm_significant']],float)
        data.update({f'pc{k}_x':x,f'pc{k}_top_lower':tl,f'pc{k}_top_median':tm,f'pc{k}_top_upper':tu,f'pc{k}_random_lower':rl,f'pc{k}_random_median':rm,f'pc{k}_random_upper':ru,f'pc{k}_significant_pct':sig,f'pc{k}_nonsignificant_pct':ns})
    return data,records
