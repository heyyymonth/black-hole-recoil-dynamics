"""Positive mixtures of finite-angle jump measures and their exact moments.

Table coordinates are radius/r_rho, speed/v_dyn and pitch in [0,1], where
pitch=2 atan2(v_t,|v_r|)/pi. Interpolation is of the entire rate measure, not
independently fitted diffusion entries. D is its raw second increment moment.
The numerical law is an approximation to the host's local encounter law;
local quadrature and state-grid refinement are separate accuracy questions.
"""
from collections import OrderedDict
from pathlib import Path
import hashlib
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for directory in (ROOT, ROOT/'analysis/host_reconstruction', ROOT/'analysis/encounters'):
    if str(directory) not in sys.path:
        sys.path.insert(0,str(directory))
from host import Host
from physics import DiscreteBath, PC_PER_KMS_MYR, scatter_impulse, velocity_ball
from analysis.encounters.run import density_scale, local_scale, elastic_floor


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def host_scales(host):
    selected = host.types == 'BH'
    core = density_scale(host,selected)
    velocity = float(np.sqrt(-core*host.psi_spline(core,1)))
    return core,velocity,PC_PER_KMS_MYR*core/velocity


def canonical_velocity(speed,pitch):
    theta = np.pi*pitch/2
    # Exact endpoint zeros make the registered symmetry cases unambiguous.
    return speed*np.array([0. if pitch==0 else np.sin(theta),0.,
                           0. if pitch==1 else np.cos(theta)])


def center_rotation(pitch):
    theta=np.pi*pitch/2
    sine=0. if pitch==0 else np.sin(theta)
    cosine=0. if pitch==1 else np.cos(theta)
    return np.array([[cosine,0.,sine],[0.,1.,0.],[-sine,0.,cosine]])


def symmetry_kind(radius,speed,pitch):
    if radius==0 and speed==0:
        return 'isotropic'
    if speed==0:
        return 'axisymmetric_even'
    if radius==0 or pitch==0:
        return 'axisymmetric'
    if pitch==1:
        return 'reflections_yz'
    return 'reflection_y'


def symmetry_moments(A,D,kind):
    """Moment average under the same exact symmetry group used by samples."""
    A,D=np.array(A,copy=True),np.array(D,copy=True)
    if kind=='isotropic':
        return np.zeros(3),np.eye(3)*np.trace(D)/3
    if kind.startswith('axisymmetric'):
        A[:2]=0.
        D=np.diag([(D[0,0]+D[1,1])/2]*2+[D[2,2]])
        if kind=='axisymmetric_even':
            A[:]=0.
        return A,D
    A[1]=0.;D[1,0]=D[0,1]=D[1,2]=D[2,1]=0.
    if kind=='reflections_yz':
        A[2]=0.;D[0,2]=D[2,0]=0.
    return A,D


def symmetry_samples(increments,kind,rng):
    increments=np.array(increments,copy=True)
    n=len(increments)
    if kind=='isotropic':
        direction=rng.normal(size=(n,3))
        direction/=np.linalg.norm(direction,axis=1)[:,None]
        return direction*np.linalg.norm(increments,axis=1)[:,None]
    if kind.startswith('axisymmetric'):
        phi=rng.uniform(0.,2*np.pi,n)
        cosine,sine=np.cos(phi),np.sin(phi)
        x,y=increments[:,0].copy(),increments[:,1].copy()
        increments[:,0]=cosine*x-sine*y
        increments[:,1]=sine*x+cosine*y
        if kind=='axisymmetric_even':
            increments[:,2]*=2*rng.integers(0,2,n)-1
    else:
        increments[:,1]*=2*rng.integers(0,2,n)-1
        if kind=='reflections_yz':
            increments[:,2]*=2*rng.integers(0,2,n)-1
    return increments


def vertex_law(host,coordinates,scales,protocol,encounter_protocol,order):
    """One positive quadrature law; physical moments are returned per Myr."""
    core,vdyn,tdyn=scales
    r,s,pitch=map(float,coordinates)
    radius,speed=r*core,s*vdyn
    base_pitch=0. if radius==0 else pitch
    velocity=canonical_velocity(speed,base_pitch)
    selected=host.types=='BH'
    bmax=protocol['bmax_over_local_scale']*local_scale(host,selected,radius,core)
    minimum=lambda w,m: elastic_floor(w,m,protocol['remnant_mass_Msun'],host.G,
                                     bmax,encounter_protocol['closest_approach'])
    relative,weights=velocity_ball(np.sqrt(2*host.psi(radius))+speed,order,
                                    radial_scale=np.sqrt(host.s2j[selected]).min())
    field=velocity-relative
    numbers=host.number_df(radius,field)[selected]*weights
    bath=DiscreteBath(field,numbers,host.mj[selected],protocol['remnant_mass_Msun'],
                      host.G,minimum,bmax)
    moment=bath.moments(velocity)
    relative_speed=np.linalg.norm(relative,axis=1)[None,:]
    lower=minimum(relative_speed,host.mj[selected,None])
    excluded=float(np.sum(numbers*relative_speed*np.pi*lower**2)/PC_PER_KMS_MYR)
    parent=float(np.sum(numbers*relative_speed*np.pi*bmax**2)/PC_PER_KMS_MYR)
    kind=symmetry_kind(radius,speed,pitch)
    A,D=symmetry_moments(moment['A'],moment['D'],kind)
    rotation=center_rotation(pitch) if radius==0 and speed>0 else np.eye(3)
    A=rotation@A;D=rotation@D@rotation.T
    exact_n=float(np.sum(host.rho_spline(radius)[selected]/host.mj[selected]))
    sigma=float(np.sqrt(np.sum((host.pr_spline(radius)[selected]
                              +host.pt_spline(radius)[selected])/host.mj[selected])/(3*exact_n)))
    answer={'rate':moment['rate'],'excluded':excluded,'parent':parent,'A':A,'D':D,
            'number_density_relative_error':abs(float(numbers.sum())/exact_n-1),
            'sigma':sigma,'kind':kind,'rotation':rotation}
    return answer,(bath,velocity,kind,rotation)


def refine_vertex(host,coordinates,scales,protocol,encounter_protocol):
    previous,history=None,[]
    limits=protocol['kernel']['local_acceptance']
    for order in protocol['kernel']['velocity_quadrature_orders']:
        answer,_=vertex_law(host,coordinates,scales,protocol,encounter_protocol,order)
        if previous is not None:
            scale=max(np.linalg.norm(answer['A']),np.trace(answer['D'])/answer['sigma'])
            differences={
                'rate':abs(answer['rate']/previous['rate']-1),
                'excluded':abs(answer['excluded']/previous['excluded']-1),
                'A':float(np.linalg.norm(answer['A']-previous['A'])/scale),
                'D':float(np.linalg.norm(answer['D']-previous['D'])/np.linalg.norm(answer['D']))}
            history.append({'order':order,**differences})
            accepted=(max(differences.values())<limits['moment_relative_refinement']
                      and answer['number_density_relative_error']<limits['number_density_relative_error']
                      and np.linalg.eigvalsh(answer['D'])[0]>=0)
            if accepted:
                break
        previous=answer
    record={'coordinates':list(map(float,coordinates)),'order':order,'accepted':bool(accepted),
            'number_density_relative_error':answer['number_density_relative_error'],
            'rate_partition_relative_error':abs((answer['rate']+answer['excluded'])/answer['parent']-1),
            'minimum_D_eigenvalue':float(np.linalg.eigvalsh(answer['D'])[0]),
            'symmetry':answer['kind'],'refinement_history':history}
    return answer,record


class VertexSampler:
    """CDF over the actual positive quadrature, followed by continuous b,phi."""
    def __init__(self,bath,velocity,kind,rotation,velocity_scale):
        relative,lower,upper=bath._state(velocity)
        rates=bath.number_weights*np.linalg.norm(relative,axis=-1)*np.pi*(upper-lower)*(upper+lower)
        positive=np.flatnonzero(rates.ravel()>0)
        if not len(positive):
            raise ValueError('A sampled vertex must have positive allowed rate')
        species,nodes=np.unravel_index(positive,rates.shape)
        self.cdf=np.cumsum(rates.ravel()[positive]);self.cdf/=self.cdf[-1];self.cdf[-1]=1.
        self.relative=relative[species,nodes]
        self.mass=bath.masses[species]
        self.lower,self.upper=lower[species,nodes],upper[species,nodes]
        self.remnant,self.G=bath.remnant_mass,bath.G
        self.kind,self.rotation,self.velocity_scale=kind,rotation,velocity_scale

    def sample(self,count,rng):
        index=np.searchsorted(self.cdf,rng.random(count),side='right')
        low,high=self.lower[index],self.upper[index]
        impact=np.sqrt(low*low+rng.random(count)*(high-low)*(high+low))
        increments=scatter_impulse(self.relative[index],self.mass[index],self.remnant,
                                   self.G,impact,rng.uniform(0.,2*np.pi,count))
        increments=symmetry_samples(increments,self.kind,rng)
        return increments@self.rotation.T/self.velocity_scale


class MatchedKernel:
    """A state-dependent positive mixture, in r_rho/v_dyn/t_dyn units.

    ``moments`` returns rate, excluded rate, A and raw D in these dimensionless
    units. ``sample_events`` returns one conditional allowed increment per
    row. Callers draw the waiting time or compound-Poisson event count.
    Querying outside either finite table axis raises; no extrapolation occurs.
    """
    def __init__(self,host_path,table_path,protocol_path,cache_size=None):
        self.protocol=json.loads(Path(protocol_path).read_text())
        self.encounter_protocol=json.loads((ROOT/self.protocol['encounter_protocol']).read_text())
        self.host=Host(dict(np.load(host_path,allow_pickle=False)))
        self.r_scale_pc,self.v_scale_kms,self.time_scale_Myr=host_scales(self.host)
        data=dict(np.load(table_path,allow_pickle=False))
        metadata=json.loads(Path(table_path).with_suffix('.json').read_text())
        if metadata['host_sha256']!=digest(host_path) or metadata['protocol_sha256']!=digest(protocol_path):
            raise ValueError('Stale kernel host or protocol')
        if metadata.get('table_sha256')!=digest(table_path):
            raise ValueError('Stale kernel table')
        if metadata['implementation_sha256']['analysis/matched/kernel.py']!=digest(Path(__file__)):
            raise ValueError('Stale kernel implementation')
        if not metadata['all_vertex_quadrature_passed']:
            raise ValueError('Unaccepted vertex quadrature cannot drive a trajectory')
        self._initialize(data)
        self._sampler_factory=self._physical_sampler
        self.cache_size=cache_size or self.protocol['kernel']['cache_vertices']
        self._cache=OrderedDict()

    @classmethod
    def from_arrays(cls,axes,rate,A,D,excluded,sampler_factory=None):
        """Dimensionless synthetic law interface for independent tests."""
        result=cls.__new__(cls)
        result._initialize({'radius':axes[0],'speed':axes[1],'pitch':axes[2],
                            'rate':rate,'A':A,'D':D,'excluded':excluded,
                            'order':np.zeros(np.shape(rate),dtype=int)})
        result._sampler_factory=sampler_factory
        result.cache_size=128;result._cache=OrderedDict()
        return result

    def _initialize(self,data):
        self.axes=[np.asarray(data[key],dtype=float) for key in ('radius','speed','pitch')]
        if any(a.ndim!=1 or len(a)<2 or not np.isfinite(a).all()
               or np.any(np.diff(a)<=0) for a in self.axes):
            raise ValueError('Increasing finite axes are required')
        if self.axes[0][0]!=0 or self.axes[1][0]!=0 or self.axes[2][0]!=0 or self.axes[2][-1]!=1:
            raise ValueError('Axes must include the origin and entire pitch interval')
        self.shape=tuple(map(len,self.axes))
        count=int(np.prod(self.shape))
        self.rate=np.asarray(data['rate'],dtype=float).reshape(count)
        self.excluded=np.asarray(data['excluded'],dtype=float).reshape(count)
        self.A=np.asarray(data['A'],dtype=float).reshape(count,3)
        self.D=np.asarray(data['D'],dtype=float).reshape(count,3,3)
        self.orders=np.asarray(data['order'],dtype=int).reshape(count)
        if (not all(np.isfinite(a).all() for a in (self.rate,self.excluded,self.A,self.D))
                or np.any(self.rate<0) or np.any(self.excluded<0)
                or not np.allclose(self.D,self.D.transpose(0,2,1),rtol=0,atol=1e-12)
                or np.min(np.linalg.eigvalsh(self.D)) < -1e-12):
            raise ValueError('Invalid positive moment table')

    def states_can_query(self,states):
        states=np.asarray(states,dtype=float)
        if states.ndim!=2 or states.shape[1]!=6:
            raise ValueError('States must be N by six')
        return (np.isfinite(states).all(axis=1)
                &(np.linalg.norm(states[:,:3],axis=1)<=self.axes[0][-1])
                &(np.linalg.norm(states[:,3:],axis=1)<=self.axes[1][-1]))

    def _coordinates(self,states):
        states=np.asarray(states,dtype=float)
        if not np.all(self.states_can_query(states)):
            raise ValueError('State outside finite kernel domain; stop rather than extrapolate')
        position,velocity=states[:,:3],states[:,3:]
        radius=np.linalg.norm(position,axis=1);speed=np.linalg.norm(velocity,axis=1)
        er=np.divide(position,radius[:,None],out=np.zeros_like(position),where=radius[:,None]>0)
        central=radius==0
        er[central]=np.divide(velocity[central],speed[central,None],
                              out=np.tile([0.,0.,1.],(central.sum(),1)),
                              where=speed[central,None]>0)
        vr=np.einsum('ni,ni->n',velocity,er)
        vt_vector=velocity-vr[:,None]*er
        vt=np.linalg.norm(vt_vector,axis=1)
        degenerate=vt<=64*np.finfo(float).eps*np.maximum(speed,1.)
        et=np.divide(vt_vector,vt[:,None],out=np.zeros_like(velocity),where=vt[:,None]>0)
        axes=np.eye(3)[np.argmin(abs(er[degenerate]),axis=1)]
        tangent=axes-np.einsum('ni,ni->n',axes,er[degenerate])[:,None]*er[degenerate]
        et[degenerate]=tangent/np.linalg.norm(tangent,axis=1)[:,None]
        eb=np.cross(er,et)
        sign=np.where(vr<0,-1.,1.)
        transform=np.stack((et,eb,sign[:,None]*er),axis=2)
        pitch=2/np.pi*np.arctan2(np.where(degenerate,0.,vt),abs(vr))
        pitch[speed==0]=0.
        return np.column_stack((radius,speed,pitch)),transform

    def _vertices(self,states):
        coordinates,transform=self._coordinates(states)
        lowers=[];fractions=[]
        for i,axis in enumerate(self.axes):
            lower=np.clip(np.searchsorted(axis,coordinates[:,i],side='right')-1,0,len(axis)-2)
            fraction=(coordinates[:,i]-axis[lower])/(axis[lower+1]-axis[lower])
            lowers.append(lower);fractions.append(fraction)
        vertex,weight=[],[]
        for dr in (0,1):
            for dv in (0,1):
                for dp in (0,1):
                    indices=(lowers[0]+dr,lowers[1]+dv,lowers[2]+dp)
                    vertex.append(np.ravel_multi_index(indices,self.shape))
                    weight.append(np.prod([f if bit else 1-f for f,bit in zip(fractions,(dr,dv,dp))],axis=0))
        return np.array(vertex).T,np.array(weight).T,transform

    def moments(self,states):
        vertex,weight,transform=self._vertices(states)
        A=np.einsum('nk,nka->na',weight,self.A[vertex])
        D=np.einsum('nk,nkab->nab',weight,self.D[vertex])
        return {'rate':np.sum(weight*self.rate[vertex],axis=1),
                'excluded':np.sum(weight*self.excluded[vertex],axis=1),
                'A':np.einsum('nij,nj->ni',transform,A),
                'D':np.einsum('nia,nab,njb->nij',transform,D,transform)}

    def _physical_sampler(self,index):
        coordinates=[self.axes[i][j] for i,j in enumerate(np.unravel_index(index,self.shape))]
        scales=(self.r_scale_pc,self.v_scale_kms,self.time_scale_Myr)
        answer,law=vertex_law(self.host,coordinates,scales,self.protocol,self.encounter_protocol,
                              int(self.orders[index]))
        for key,factor in [('rate',self.time_scale_Myr),('excluded',self.time_scale_Myr),
                           ('A',self.time_scale_Myr/self.v_scale_kms),
                           ('D',self.time_scale_Myr/self.v_scale_kms**2)]:
            if not np.allclose(np.asarray(answer[key])*factor,getattr(self,key)[index],rtol=2e-11,atol=2e-13):
                raise ValueError('Rebuilt vertex law differs from its recorded moments')
        return VertexSampler(*law,self.v_scale_kms)

    def _sampler(self,index):
        index=int(index)
        if index in self._cache:
            self._cache.move_to_end(index)
        else:
            if self._sampler_factory is None:
                raise ValueError('No conditional sampler supplied')
            self._cache[index]=self._sampler_factory(index)
            if len(self._cache)>self.cache_size:
                self._cache.popitem(last=False)
        return self._cache[index]

    def sample_events(self,states,rng):
        """One independent allowed event per supplied state, in global axes."""
        vertex,weight,transform=self._vertices(states)
        flux=weight*self.rate[vertex]
        totals=flux.sum(axis=1)
        if np.any(totals<=0):
            raise ValueError('Cannot sample events at zero-rate states')
        cumulative=np.cumsum(flux/totals[:,None],axis=1);cumulative[:,-1]=1.
        choice=np.sum(rng.random(len(vertex))[:,None]>=cumulative,axis=1)
        selected=vertex[np.arange(len(vertex)),choice]
        increments=np.empty((len(vertex),3))
        for index in np.unique(selected):
            rows=np.flatnonzero(selected==index)
            increments[rows]=self._sampler(index).sample(len(rows),rng)
        return np.einsum('nij,nj->ni',transform,increments)
