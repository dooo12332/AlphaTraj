#%%
import os
import sys
import argparse
import alphaspace2 as al
from alphaspace2 import Snapshot as ss
from typing import List,Union,Dict,Tuple
import mdtraj
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import fcluster, linkage
from configparser import ConfigParser
from tqdm import tqdm, trange
#%% process class

class SnapShots:
    def __init__(self,rec:mdtraj.Trajectory) -> None:
        self.rec:mdtraj.Trajectory=rec
        #al.annotateVinaAtomTypes(pdbqt=rec_pdbqt, receptor=self.rec)
        self.sslist:List[al.Snapshot.Snapshot]=[]

        self.boxCenterAtomId:List[List[int]]=[]
        self.boxlen:List[List[float]]=[]
        
        self.distance_cutoff=3.5

    def __getitem__(self,id:int):
        return self.sslist[id]

    def SetBox(self,center_atom_id:List[List[int]],length:List[List[float]]):
        self.boxCenterAtomId=center_atom_id
        self.boxlen=length

    def _CalcBoxEdge(self,frame:int)->np.ndarray:
        xyz=np.squeeze(self.rec[frame].xyz*10,0)
        box_num=len(self.boxCenterAtomId)
        out=np.zeros((box_num,6))
        for i in range(box_num):
            txyz=np.average(xyz[np.array(self.boxCenterAtomId[i])-1],axis=0)
            out[i,0]=txyz[0]-self.boxlen[i][0]/2.
            out[i,1]=txyz[0]+self.boxlen[i][0]/2.
            out[i,2]=txyz[1]-self.boxlen[i][1]/2.
            out[i,3]=txyz[1]+self.boxlen[i][1]/2.
            out[i,4]=txyz[2]-self.boxlen[i][2]/2.
            out[i,5]=txyz[2]+self.boxlen[i][2]/2.
        return out

    def Process(self,start:int=0,stop:int=-1,offset:int=1,rec_mask:np.ndarray=np.array([]),lig_mask:np.ndarray=np.array([])):
        if start<0:
            start=0
        elif  start>self.rec.n_frames:
            start=self.rec.n_frames

        if stop<start:
            stop=self.rec.n_frames
        elif stop>self.rec.n_frames:
            stop=self.rec.n_frames
        
        print(f'start generate snapshots list({start}:{stop}:{offset}). This may take a long time...')
        for i in tqdm(range(start,stop,offset),desc='Processing',unit='frame',unit_scale=True):
            self.sslist.append(al.Snapshot())
            if len(self.boxCenterAtomId)!=0:
                #print(f'boxedge={self._CalcBoxEdge(i)}')
                self.sslist[-1].SetBox(self._CalcBoxEdge(i))
            self.sslist[-1].run(receptor=self.rec,frame=i,rec_mask=rec_mask,lig_mask=lig_mask)

    def GetPockComposition(self,frame:int,pid:int)->List[int]:
        return self.sslist[frame]._pocket_alpha_index_list[pid]

    def GetPockNum(self,frame:int)->int:
        return len(self.sslist[frame]._pocket_space)
    
class PocketsInfo:
    def __init__(self) -> None:
        self._pockets_cluter_id:np.ndarray=np.array([])
        self._pockets_frequency:np.ndarray
        self._pockets_score:np.ndarray
        self._pocket_cluster_splited:List[np.ndarray]=[]
        self._pocket_volume_v_time:np.ndarray#majtrix:pocket num * frame
        self._pocket_volume_v_time_normed:np.ndarray
        self._pocket_volume_v_time_sorted:np.ndarray#small->large
        self._pocket_volume_range:np.ndarray
        self._pocket_nonpolar_v_time:np.ndarray
        self._pocket_occupancy_v_time:np.ndarray
        self._pocket_occupancy_score:np.ndarray
        self._pocket_composition:Dict[str, list]={}#main_name:[main_name,main_score,minor_dic,frame_list],dic={sub_name:[sub_name,sub_score,[exists frames]]}
        self._pocket_composition_v_time:np.ndarray
        self._pocket_main_group_v_time:List[List[int]]=[]
        self._pocket_sub_group_v_time:List[List[int]]=[]
        self.fluc_range:np.ndarray

    def _CalcPocketScore(self)->None:
        result=[]
        for i in range(self.GetPocketNum()):
            #result.append(self.GetLifeTime(i)*self.GetMedianVol(i)*self.GetVolumeFluctuations(i,80)[3]/self.GetFrameNum()/self.GetMaxVol(i)*100)
            result.append(self.GetLifeTime(i)*self.GetMedianVol(i)*(1-self.GetVolumeFluctuations(i,80)[3])*(1-np.sqrt(self.GetMSENPR(i,80)))/self.GetFrameNum())
        self._pockets_score=np.array(result)

    def _CalaPocketKey(self,frame:int,main_score_cutoff:float=20.0,sub_score_cutoff:float=10.0)->Tuple[List[int],float,List[int],float]:
        main_key=[]
        sub_key=[]
        main_score=0
        sub_score=0
        tmp_comp=self.GetPocketComposition(frame)
        for j in tmp_comp:
            if self.GetScore(j)>main_score_cutoff:
                main_key.append([j,self.GetScore(j)])
                main_score+=self.GetScore(j)
            elif self.GetScore(j)>sub_score_cutoff:
                sub_key.append([j,self.GetScore(j)])
                sub_score+=self.GetScore(j)
        main_key=sorted(main_key,key=lambda x:x[1],reverse=True)
        sub_key=sorted(sub_key,key=lambda x:x[1],reverse=True)
        return ([i[0] for i in main_key],main_score,[i[0] for i in sub_key],sub_score)

    def _key2str(self,key:List[int])->str:
        if len(key)==0:
            return 'none'
        ukey=np.unique(key)
        out=''
        for i in ukey:
            out+=f'{i:03d}-'
        return out[:-1]

    def _PocketGrouppAnalysis(self,main_score_cutoff:float=20.0,sub_score_cutoff:float=10.0):
        for i in range(self.GetFrameNum()):
            #print(f'fram[{i}]    ')
            main_key,main_score,sub_key,sub_score=self._CalaPocketKey(i,main_score_cutoff,sub_score_cutoff)
            main_name=self._key2str(main_key)
            sub_name=self._key2str(sub_key)
            self._pocket_main_group_v_time.append(main_key)
            self._pocket_sub_group_v_time.append(sub_key)
            #print(f'main={main_name}, sub={sub_name}\n')
            if main_name not in self._pocket_composition.keys():
                self._pocket_composition[main_name]=[main_name,main_score,{sub_name:[sub_name,sub_score,[i]]},[i]]
            elif sub_name not in self._pocket_composition[main_name][2].keys():
                self._pocket_composition[main_name][3].append(i)
                self._pocket_composition[main_name][2][sub_name]=[sub_name,sub_score,[i]]
            else:
                self._pocket_composition[main_name][3].append(i)
                self._pocket_composition[main_name][2][sub_name][2].append(i)

    def _GetBestConf(self,frame_list:List[int],reserve_num:int,percentile:int=80,life_time_cutoff:int=10,is_use_scoring:bool=True):
        score=self._pockets_score.copy()
        if is_use_scoring==False:
            score=score*0+1
        out=[]#[[frame index,total vol/score,main group,minor group]]
        for i in range(reserve_num):
            out.append([-1,0])
        for i in frame_list:
            tsumscore=0
            tsumvol=0
            for j in range(self.GetPocketNum()):
                tvol=self._pocket_volume_v_time[j,i]
                if tvol>0 and self.GetLifeTime(j)>life_time_cutoff:
                    if tvol > self.GetPCTLVol(j,percentile):
                        tvol=self.GetPCTLVol(j,percentile)
                    tsumscore+=tvol*score[j]
                    tsumvol+=tvol
            out.append([i,tsumscore,tsumvol,self.GetMainName(i),self.GetMinorName(i)])
            # for k in range(5):
            #     if tsumscore>out[k][1]:
            #         out.insert(k,[i,tsumscore,tsumvol,self.GetMainName(i),self.GetSubName(i)])
            #         break
        out=sorted(out,key=lambda x:x[1],reverse=True)
        return out[:reserve_num]

    def GetFrameNum(self)->int:
        return self._pocket_volume_v_time.shape[1]

    def GetPocketNum(self)->int:
        return self._pockets_frequency.shape[0]

    def GetPocketComposition(self,frame:int)->np.ndarray:
        return np.unique(self._pocket_cluster_splited[frame])

    def GetMainGroup(self,frame:int)->List[int]:
        return self._pocket_main_group_v_time[frame]
    
    def GetMainName(self,frame:int)->str:
        return self._key2str(self.GetMainGroup(frame))

    def GetMinorGroup(self,frame:int)->List[int]:
        return self._pocket_sub_group_v_time[frame]
    
    def GetMinorName(self,frame:int)->str:
        return self._key2str(self.GetMinorGroup(frame))
    
    def GetGroupName(self,frame:int)->Tuple[List[int],List[int]]:
        return (self.GetMainGroup(frame),self.GetMinorGroup(frame))

    def GetLifeTime(self,pid:int)->float:
        return self._pockets_frequency[pid]

    def GetAllLifeTime(self)->np.ndarray:
        return self._pockets_frequency.copy()

    def GetVol(self,pid:int,frame:int)->float:
        return self._pocket_volume_v_time[pid,frame]
    
    def HistogramVol(self,pid:int,bin_num:int)->Tuple[np.ndarray,np.ndarray]:
        tvol=self._pocket_volume_v_time_sorted[pid]
        tvol=tvol[tvol>0]
        return np.histogram(tvol,bin_num)

    def GetMaxVol(self,pid:int)->float:
        return self._pocket_volume_v_time_sorted[pid,-1]
    
    def GetMinVol(self,pid:int)->float:
        tvol=self._pocket_volume_v_time_sorted[pid]
        tvol=tvol[tvol>0]
        return tvol[0]

    def GetPCTLVol(self,pid:int,percentile:int)->float:
        tvol=self._pocket_volume_v_time_sorted[pid]
        tvol=tvol[tvol>0]
        return np.percentile(tvol,percentile)

    def GetMeanVol(self,pid:int):
        tvol=self._pocket_volume_v_time_sorted[pid]
        tvol=tvol[tvol>0]
        return np.mean(tvol)
    
    def GetMedianVol(self,pid:int)->float:
        tvol=self._pocket_volume_v_time_sorted[pid]
        tvol=tvol[tvol>0]
        return np.median(tvol)

    def GetMSEVol(self,pid:int)->float:
        tvol=self._pocket_volume_v_time_sorted[pid]
        tvol=tvol[tvol>0]
        return np.mean(np.power((tvol-np.mean(tvol)),2))

    def GetVolumeFluctuations(self,pid:int,percentage:int)->Tuple[float,float,float,float]:
        tvol=self._pocket_volume_v_time_sorted[pid]
        tvol=tvol[tvol>0]
        if tvol.shape[0]==0:
            return (0,0,0,0)
        if tvol.shape[0]==1:
            return (tvol[0],tvol[0],0,0)
        if tvol.shape[0]<10:
            tmin=tvol[0]
            tmax=tvol[-1]
            #return (tvol[0],tvol[0],0,0)
        else:
            tmin=np.percentile(tvol,int((100-percentage)/2))
            tmax=np.percentile(tvol,int((50+percentage/2)))
        return (tmin,tmax,tmax-tmin,(tmax-tmin)/(tvol[-1]-tvol[0]+1e-5))

    def GetVolVTime(self,pid:int)->np.ndarray:
        return self._pocket_volume_v_time[pid].copy()

    def GetTotalVolVTime(self,score_cutoff=10)->np.ndarray:
        return np.sum(self._pocket_volume_v_time[self._pockets_score>score_cutoff],axis=0)

    def GetNoZeroVolVTime(self,pid:int)->np.ndarray:
        tvol=self._pocket_volume_v_time_sorted[pid]
        return tvol[tvol>0]

    def GetNPR(self,pid:int,frame:int)->float:
        if self._pocket_nonpolar_v_time[pid,frame]==0:
            return float('nan')
        return self._pocket_nonpolar_v_time[pid,frame]

    def GetMaxNPR(self,pid:int)->float:
        tnpr=self._pocket_nonpolar_v_time[pid]
        tnpr=tnpr[tnpr>0]
        return np.max(tnpr)

    def GetMinNPR(self,pid:int)->float:
        tnpr=self._pocket_nonpolar_v_time[pid]
        tnpr=tnpr[tnpr>0]
        return np.min(tnpr)

    def GetMeanNPR(self,pid:int,percentile:int=100)->float:
        sort_id=np.argsort(self._pocket_volume_v_time[pid])
        tnpr=self._pocket_nonpolar_v_time[pid,sort_id]
        tnpr=tnpr[self._pocket_volume_v_time_sorted[pid]>0]
        l=tnpr.shape[0]
        if l<10:
            return np.mean(tnpr)
        return np.mean(tnpr[int((100-percentile)/200*l):int(np.round((0.5+percentile/200)*l))])
        
    def GetMSENPR(self,pid:int,percentile:int=100)->float:
        sort_id=np.argsort(self._pocket_volume_v_time[pid])
        tnpr=self._pocket_nonpolar_v_time[pid,sort_id]
        tnpr=tnpr[self._pocket_volume_v_time_sorted[pid]>0]
        l=tnpr.shape[0]
        if l<10:
            return np.mean(np.power((tnpr-np.mean(tnpr)),2))
        tnpr=tnpr[int((100-percentile)/200*l):int(np.round((0.5+percentile/200)*l))]
        return np.mean(np.power((tnpr-np.mean(tnpr)),2))

    def GetNPRVTime(self,pid:int)->np.ndarray:#0==nan
        return self._pocket_nonpolar_v_time[pid].copy()

    def GetNoZeroNPRVTime(self,pid:int)->np.ndarray:
        tnpr=self._pocket_nonpolar_v_time[pid]
        tnpr=tnpr[tnpr>0]
        return tnpr.copy()

    def GetNPRFluctuations(self,pid:int,percentage:int)->Tuple[float,float,float,float]:
        tnpr=self._pocket_nonpolar_v_time[pid]
        tnpr=tnpr[tnpr>0]
        if tnpr.shape[0]==1:
            return (tnpr[0],tnpr[0],0,0)
        tnpr=np.sort(tnpr)
        tmin=np.percentile(tnpr,int((100-percentage)/2))
        tmax=np.percentile(tnpr,int((50+percentage/2)))
        return (tmin,tmax,tmax-tmin,(tmax-tmin)/(tnpr[-1]-tnpr[0]+1e-5))

    def GetOccupancy(self,pid:int,frame:int)->float:
        return self._pocket_occupancy_v_time[pid,frame]

    def GetOccupancyRatio(self,pid:int,frame:int)->float:
        return self.GetOccupancy(pid,frame)/(self.GetVol(pid,frame)+1e-5)

    def GetMeanOccRatio(self,pid:int,percentile:int=80)->float:
        occ_ratio=self._pocket_occupancy_v_time[pid]/(self._pocket_volume_v_time[pid]+1e-5)
        sort_id=np.argsort(self._pocket_volume_v_time[pid])
        occ_ratio_s=occ_ratio[sort_id]
        occ_ratio_s=occ_ratio_s[self._pocket_volume_v_time_sorted[pid]>0]
        l=occ_ratio_s.shape[0]
        if l<10:
            return np.mean(occ_ratio_s)
        return np.mean(occ_ratio_s[int((100-percentile)/200*l):int(np.round((0.5+percentile/200)*l))])

    def GetOccupancyVTime(self,pid:int)->np.ndarray:
        return self._pocket_occupancy_v_time[pid].copy()

    def GetScore(self,pid:int)->float:
        return self._pockets_score[pid]
    
    def GetAllScore(self)->np.ndarray:
        return self._pockets_score.copy()

    def SyncBTWPockets(self)->np.ndarray:
        pnum=self.GetPocketNum()
        out=np.eye(pnum)
        for i in range(pnum):
            for j in range(i+1,pnum):
                pivt=self.GetVolVTime(i)
                pivt[pivt>0]=2
                povt=self.GetVolVTime(j)
                povt[povt>0]=2
                total_frame=pivt+povt
                total_frame=total_frame[total_frame>0].shape[0]+1e-5
                intersection=pivt*povt
                intersection=intersection[intersection==4].shape[0]
                out[i,j]=out[j,i]=intersection/total_frame
        return out

    def PocketAcorrAnalysis(self):
        pnum=self.GetPocketNum()
        out=np.eye(pnum)
        for i in range(pnum):
            for j in range(i+1,pnum):
                pivt=(self.GetVolVTime(i)[1::]-self.GetVolVTime(i)[:-1:])/(self.GetMaxVol(i)+1e-5)
                povt=(self.GetVolVTime(j)[1::]-self.GetVolVTime(j)[:-1:])/(self.GetMaxVol(j)+1e-5)
                out[i,j]=out[j,i]=np.sum(pivt*povt)/(np.linalg.norm(pivt)*np.linalg.norm(povt)+1e-5)
        return out

    def GetBestConf(self,reserve_num:int,percentile:int=80,life_time_cutoff:int=10,is_use_scoring:bool=True):
        frame_list=list(range(self.GetFrameNum()))
        return self._GetBestConf(frame_list,reserve_num,percentile,life_time_cutoff,is_use_scoring)

    def GetMainGroupName(self)->List[str]:
        return list(self._pocket_composition.keys())

    def GetMainGroupNum(self)->int:
        return len(self._pocket_composition)

    def GetMainGroupFrames(self,main_name)->int:
        return len(self._pocket_composition[main_name][-1])

    def GetMinorGroupName(self,main_name)->List[str]:
        return list(self._pocket_composition[main_name][2].keys())

    def GetMinorGroupFrames(self,main_name,minor_name)->int:
        return len(self._pocket_composition[main_name][2][minor_name][-1])

    def GetMinorGroupNum(self,main_name:str)->int:
        out=self._pocket_composition.get(main_name,None)
        if out==None:
            return 0
        return len(out[2])

    def GetMainGroupVTime(self)->List[str]:
        return list(map(self._key2str,self._pocket_main_group_v_time))

    def GetAllPocketGroupVTime(self)->Tuple[list,np.ndarray]:
        name=[]
        size=0
        count=0
        for k,v in self._pocket_composition.items():
            size+=len(v[2])
        self._pocket_composition_v_time=np.zeros((size,self.GetFrameNum()))-1
        main_list=list(self._pocket_composition.values())
        main_list=sorted(main_list,key=lambda x:x[1],reverse=True)
        for i in main_list:
            sub_list=list(i[2].values())
            if len(sub_list)>1:
                sub_list=sorted(sub_list,key=lambda x:x[1],reverse=True)
            for j in sub_list:
                name.append(f'{i[0]},{j[0]}')
                self._pocket_composition_v_time[count,j[2]]=size-count
                count+=1
        return (name,self._pocket_composition_v_time.copy())

    def GetMainBestConf(self,main_name:str,reserve_num:int,percentile:int=80,life_time_cutoff:int=10,is_use_scoring:bool=True):
        return self._GetBestConf(self._pocket_composition[main_name][3],reserve_num,percentile,life_time_cutoff,is_use_scoring)

    def GetGroupBestConf(self,main_name:str,sub_name:str,reserve_num:int,percentile:int=80,life_time_cutoff:int=10,is_use_scoring:bool=True):
        return self._GetBestConf(self._pocket_composition[main_name][2][sub_name][2],reserve_num,percentile,life_time_cutoff,is_use_scoring)

    def MainPocketTransProbAnalysis(self,frame_cutoff=0)->Tuple[List[str],np.ndarray]:
        main_list=list(self._pocket_composition.values())
        main_list=sorted(main_list,key=lambda x:x[1],reverse=True)#sort by score
        main_frames_num_list=np.array([len(i[-1]) for i in main_list])
        main_id=np.argwhere(main_frames_num_list>frame_cutoff)
        main_list=[main_list[i] for i in main_id.flatten()]
        mname_list=[i[0] for i in main_list]
        out=np.zeros((len(main_list),len(main_list)))
        m_v_t=self.GetMainGroupVTime()
        for i in range(len(m_v_t)-1):
            if m_v_t[i]!=m_v_t[i+1]:
                out[mname_list.index(m_v_t[i]),mname_list.index(m_v_t[i+1])]+=1
        return (mname_list,out)

class PocketsAnalysis:
    def __init__(self,rec:mdtraj.Trajectory,rec_mask:List[int]=[],lig_mask:List[int]=[]) -> None:
        self.rec:mdtraj.Trajectory=rec
        self.lig_mask:np.ndarray=np.array(lig_mask)
        self.rec_mask:np.ndarray=np.array(rec_mask)
        if self.lig_mask.shape[0]!=0:
            self.lig_mask-=1
        if self.rec_mask.shape[0]==0:
            self.rec_mask:np.ndarray=np.arange(0,self.rec.xyz.shape[1])
            self.rec_mask=np.setdiff1d(self.rec_mask,self.lig_mask,True)
        else:
            self.rec_mask-=1
        
        rec=rec.superpose(rec,0)
        self.snap_shots=SnapShots(rec)
        self.pocketsinfo=PocketsInfo()
        self.distance_cutoff:float=3.0

    def size(self)->int:
        '''
        return frame number
        '''
        return len(self.snap_shots.sslist)
    
    def _ClusterPockets(self):
        total_coords=[]
        #Extract pocket coordinates
        for ss in self.snap_shots.sslist:
            total_coords.extend(ss._pocket_xyz)
        total_coords=np.array(total_coords)
        #pocket clustering
        self.pocketsinfo._pockets_cluter_id=fcluster(linkage(total_coords, method='average'), self.distance_cutoff,criterion='distance')
        self.pocketsinfo._pockets_cluter_id-=1

    def _PostProcessing(self)->None:
        pnum_per_frame=[]
        #Extract pocket frame num
        for ss in self.snap_shots.sslist:
            pnum_per_frame.append(len(ss._pocket_space))
        self.pocketsinfo._pockets_frequency=np.bincount(self.pocketsinfo._pockets_cluter_id)
        #generate pocketvolume_vs_time array & Delete frequency of multiple statistics
        cumsum_pnum_per_frame=np.cumsum(pnum_per_frame)
        self.pocketsinfo._pocket_cluster_splited=np.split(self.pocketsinfo._pockets_cluter_id,cumsum_pnum_per_frame)
        self.pocketsinfo._pocket_volume_v_time=np.zeros((np.max(self.pocketsinfo._pockets_cluter_id)+1,len(self.snap_shots.sslist)))
        self.pocketsinfo._pocket_nonpolar_v_time=np.zeros((np.max(self.pocketsinfo._pockets_cluter_id)+1,len(self.snap_shots.sslist)))
        self.pocketsinfo._pocket_occupancy_v_time=np.zeros((np.max(self.pocketsinfo._pockets_cluter_id)+1,len(self.snap_shots.sslist)))
        for i in range(len(self.snap_shots.sslist)):
            tpids=self.pocketsinfo._pocket_cluster_splited[i]#tmp pockets index
            tpids_count=np.bincount(tpids)
            for j in np.argwhere(tpids_count>1):#Delete frequency of multiple statistics
                self.pocketsinfo._pockets_frequency[j]-=1
            #print(f'frame[{i}] overlap {self.pocketsinfo._pockets_frequency[np.argwhere(tpids_count>1)]}')
            for j in range(tpids.shape[0]):
                self.pocketsinfo._pocket_volume_v_time[tpids[j],i]+=self.snap_shots.sslist[i]._pocket_space[j]
                self.pocketsinfo._pocket_nonpolar_v_time[tpids[j],i]+=np.mean(self.snap_shots.sslist[i]._alpha_space_nonpolar_ratio[self.snap_shots.sslist[i]._pocket_alpha_index_list[j]])+1e-5
        self.pocketsinfo._pocket_volume_v_time_sorted=np.sort(self.pocketsinfo._pocket_volume_v_time)
        if  self.rec_mask.shape[0]!=0 and self.lig_mask.shape[0]!=0:
            for i in range(len(self.snap_shots.sslist)):
                tpids=self.pocketsinfo._pocket_cluster_splited[i]
                for j in range(tpids.shape[0]):
                    taids=self.snap_shots[i]._pocket_alpha_index_list[j]
                    taspace=self.snap_shots[i]._alpha_space[taids]
                    tiscontact=self.snap_shots[i]._alpha_contact[taids]
                    self.pocketsinfo._pocket_occupancy_v_time[tpids[j],i]+=np.sum(taspace[tiscontact])
        
    def SetBox(self,center_atom_id:List[List[int]],length:List[List[float]]):
        if len(center_atom_id)!=len(length):
            print('Input parameter dimension mismatch!')
        else:
            is_set=True
            for i in length:
                if len(i)!=3:
                    print('Dimension mismatch of parameter length')
                    break
            if is_set:
                self.snap_shots.SetBox(center_atom_id,length)

    def SetDistCutoff(self,distance:float):
        self.distance_cutoff=distance

    def Analysis(self,start:int=0,end:int=-1,offset:int=1):
        print('processing snap shots....     ')
        self.offset=offset
        self.snap_shots.Process(start,end,offset,self.rec_mask,self.lig_mask)
        print('Clustering pockets...    ',end='')
        self._ClusterPockets()
        print(' done.\nAnalysing pockets...    ',end='')
        self._PostProcessing()
        print(' done.\nScoring pockets...   ',end='')
        self.pocketsinfo._CalcPocketScore()
        self.pocketsinfo._PocketGrouppAnalysis()
        print(' done.\n')

    def WritePockets(self,frame:int,out_file:str)->None:
        resi=1
        serial=1
        piday=self.pocketsinfo._pocket_cluster_splited[frame]
        pxyz=self.snap_shots.sslist[frame]._pocket_xyz
        axyz=self.snap_shots.sslist[frame]._alpha_xyz
        #print(f'piday={piday}')
        with open(out_file,'w') as f:
            tline=''
            for indx in range(piday.shape[0]):
                #print(f'indx={indx},pockcomp={self.snap_shots.GetPockComposition(frame,indx)}')
                #serial[]atom_name[]res_name[][2]res_seq[][][][4]xyz-occupancy-tempfactor[][][][][][][][][][][11]element-charge
                tline=f"{'ATOM':<6s}{serial:>5d} {'PCC':^4s} {'APL':>3s}  {resi:>4d}    {pxyz[indx][0]:8.3f}{pxyz[indx][1]:8.3f}{pxyz[indx][2]:8.3f}{1.0:6.2f}{piday[indx]:6.2f}          {'PC':2s} 0\n"
                f.writelines(tline)
                serial+=1
                for aid in self.snap_shots.GetPockComposition(frame,indx):
                    #print(f'aid={aid}')
                    if self.snap_shots.sslist[frame]._alpha_contact[aid]==0:
                        tline=f"{'ATOM':<6s}{serial:>5d} {'AAU':^4s} {'APL':>3s}  {resi:>4d}    {axyz[aid,0]:8.3f}{axyz[aid,1]:8.3f}{axyz[aid,2]:8.3f}{1.0:6.2f}{piday[indx]:6.2f}          {'AC':2s} 0\n"
                    else:
                        tline=f"{'ATOM':<6s}{serial:>5d} {'AAO':^4s} {'APL':>3s}  {resi:>4d}    {axyz[aid,0]:8.3f}{axyz[aid,1]:8.3f}{axyz[aid,2]:8.3f}{1.0:6.2f}{piday[indx]:6.2f}          {'AC':2s} 0\n"
                    f.writelines(tline)
                    serial+=1
                resi+=1

class ModelGroup:
    def __init__(self,**kw) -> None:
        self.pa_list:List[PocketsAnalysis]=[]
        self.distance_cutoff=kw.get('disCutoff',3.0)

    def __getitem__(self,id:int)->PocketsAnalysis:
        return self.pa_list[id]

    def SetDistCutoff(self,distance:float):
        self.distance_cutoff=distance
        for pa in self.pa_list:
            pa.SetDistCutoff(distance)

    def size(self)->int:
        return len(self.pa_list)

    def AddModel(self,rec:mdtraj.Trajectory,rec_mask:List[int],lig_mask:List[int])->None:
        pa=PocketsAnalysis(rec,rec_mask,lig_mask)
        pa.SetDistCutoff(self.distance_cutoff)
        self.pa_list.append(pa)

    def AddPA(self,pa:PocketsAnalysis)->None:
        pa.SetDistCutoff(self.distance_cutoff)
        self.pa_list.append(pa)

    def GetModel(self,sid:int)->PocketsAnalysis:
        return self.pa_list[sid]

    def AlignModel(self,rec_mask:List[List[int]])->None:
        '''
        rec_mask:[ [system1 atom indexs],   [system2 atom indexs] .....   ]
        if len(rec_mask)==1 : All systems will apply the same mask
        '''
        for i in range(1,self.size()):
            self.pa_list[i].snap_shots.rec.superpose(self.pa_list[0].snap_shots.rec,frame=0,ref_atom_indices=rec_mask[0],atom_indices=rec_mask[i])

    def _ClusterPockets(self)->None:
        total_coords=[]
        #Extract pocket coordinates
        for pa in self.pa_list:
            for ss in pa.snap_shots.sslist:
                total_coords.extend(ss._pocket_xyz)
        total_coords=np.array(total_coords)
        #pocket clustering
        pockets_cluter_id=fcluster(linkage(total_coords, method='average'), self.distance_cutoff,criterion='distance')
        pockets_cluter_id-=1
        system_frame_num=[pa.size() for pa in self.pa_list]
        pockets_cluter_id_list=np.split(pockets_cluter_id,np.cumsum(system_frame_num))
        for i in range(self.size()):
            #print(f'len={len(pockets_cluter_id_list[i])}')
            self.pa_list[i].pocketsinfo._pockets_cluter_id=pockets_cluter_id_list[i]

    def Analysis(self,step:List[List[int]])->None:
        for i in range(self.size()):
            print(f'processing snap shots of MODEL{i}....     ')
            self.pa_list[i].offset=step[i][2]
            self.pa_list[i].snap_shots.Process(step[i][0],step[i][1],step[i][2],self.pa_list[i].rec_mask,self.pa_list[i].lig_mask)
            print(' done.\n')
        print('Clustering pockets...    ',end='')
        self._ClusterPockets()
        print(' done.')
        for i in range(self.size()):
            print(f'Analysing pockets of MODEL{i}...    ',end='')
            self.pa_list[i]._PostProcessing()
            print(f' done.\nScoring pockets of MODEL{i}...   ',end='')
            self.pa_list[i].pocketsinfo._CalcPocketScore()
            self.pa_list[i].pocketsinfo._PocketGrouppAnalysis()
            print(' done.\n')

#%% helper class
class PAHelper:
    @staticmethod
    def _WriteBestFrames(pa:PocketsAnalysis,out_dir:str,result:list,**kw)->None:
        frame_shift=kw.get('frame_shift',0)
    
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        pdb_dir=os.path.join(out_dir,'pdb')
        if not os.path.exists(pdb_dir):
            os.makedirs(pdb_dir)

        if pa.lig_mask.shape[0]!=0:
            mean_occu=[]
            for i in range(pa.pocketsinfo.GetPocketNum()):
                mean_occu.append(pa.pocketsinfo.GetMeanOccRatio(i,80))
        else:
            mean_occu=np.zeros(pa.pocketsinfo.GetPocketNum())

        all_p_rank=[[i,pa.pocketsinfo.GetScore(i)] for i in range(pa.pocketsinfo.GetPocketNum()) ]
        all_p_rank=np.array(sorted(all_p_rank,key=lambda x:x[1],reverse=True))[:,0].transpose()

        out=result
        info=[]
        for i in range(len(out)):
            pcom=pa.pocketsinfo.GetPocketComposition(out[i][0])
            tmp=[[pcom[j],pa.pocketsinfo.GetScore(pcom[j]),pa.pocketsinfo.GetVol(pcom[j],out[i][0])]for j in range(len(pcom))]
            tmp=sorted(tmp,key=lambda x:x[1],reverse=True)
            info.append(tmp)
            pa.WritePockets(out[i][0],os.path.join(pdb_dir,f'apocket-r{i:03d}.pdb'))
            pa.rec[frame_shift+out[i][0]].save(os.path.join(pdb_dir,f'protein-r{i:03d}-f{frame_shift+out[i][0]*pa.offset:05d}.pdb'))
        with open(os.path.join(out_dir,'out.dat'),'w') as f:
            for i in range(len(out)):
                if out[i][0]==-1:
                    continue
                f.writelines(f'rank={i+1} frame={out[i][0]*pa.offset+frame_shift+1:6d} pock_score={out[i][1]:10.3f} pock_space={out[i][2]} main_name={out[i][3]} minor_name={out[i][4]}\n{"id":^3s} {"rank":^4s} {"socre":^6s} {"space":^7s} {"mnpr%80":^8s} {"occu":^6s} {"m_occu":^6s}\n')
                for j in info[i]:
                    f.writelines(f'{j[0]:3d} {np.argwhere(all_p_rank==j[0])[0,0]+1:4d} {j[1]:6.2f} {j[2]:7.2f} {pa.pocketsinfo.GetMeanNPR(j[0],80)*100:8.2f} {pa.pocketsinfo.GetOccupancyRatio(j[0],out[i][0]):6.4f} {mean_occu[j[0]]:6.4f}\n')
                f.writelines('\n\n')

        np.savetxt(os.path.join(out_dir,'voltotal_v_time.dat'),pa.pocketsinfo.GetTotalVolVTime())
        
    @staticmethod
    def WriteBestFrames(pa:PocketsAnalysis,out_dir:str,**kw)->None:
        '''
        reserve_num(type:int default:10)
        percentile(type:int range:0-100 default:80)
        life_time_cutoff(type:int default:10)
        is_use_scoring(type:bool default:True)
        frame_shift(type:int default:0)
        '''
        reserve_num=kw.get('reserve_num',10)
        percentile=kw.get('percentile',80)
        life_time_cutoff=kw.get('life_time_cutoff',10)
        is_use_scoring=kw.get('is_use_scoring',True)
        frame_shift=kw.get('frame_shift',0)
        #score_cutoff=kw.get('score_cutoff',10)

        out=pa.pocketsinfo.GetBestConf(reserve_num=reserve_num,percentile=percentile,life_time_cutoff=life_time_cutoff,is_use_scoring=is_use_scoring)
        PAHelper._WriteBestFrames(pa,out_dir,out,frame_shift=frame_shift)

    @staticmethod
    def WriteMainBestFrames(pa:PocketsAnalysis,out_dir:str,**kw)->None:
        '''
        reserve_num(type:int default:10)
        percentile(type:int range:0-100 default:80)
        life_time_cutoff(type:int default:10)
        is_use_scoring(type:bool default:True)
        frame_shift(type:int default:0)
        '''
        reserve_num=kw.get('reserve_num',1)
        percentile=kw.get('percentile',80)
        life_time_cutoff=kw.get('life_time_cutoff',10)
        is_use_scoring=kw.get('is_use_scoring',True)
        frame_shift=kw.get('frame_shift',0)

        main_list=list(pa.pocketsinfo._pocket_composition.values())
        main_list=sorted(main_list,key=lambda x:x[1],reverse=True)
        for i in main_list:
            m_out_dir=os.path.join(out_dir,f'main_{i[0]}_s{np.round(i[1])}')
            data=pa.pocketsinfo.GetMainBestConf(i[0],reserve_num=reserve_num,percentile=percentile,life_time_cutoff=life_time_cutoff,is_use_scoring=is_use_scoring)
            PAHelper._WriteBestFrames(pa,m_out_dir,data,frame_shift=frame_shift)

    @staticmethod
    def WritePocketAccor(pa:PocketsAnalysis,out_dir:str,**kw)->None:
        score_cutoff=kw.get('score_cutoff',10.0)
        sync_name=kw.get('coexist_name','Pocket_coexistence rate.dat')
        accor_name=kw.get('accor_name','Pocket_correlation.dat')
        id_name=kw.get('id_name','Pocket_sync_corr_id.dat')

        m=pa.pocketsinfo.SyncBTWPockets()
        m=m[pa.pocketsinfo.GetAllScore()>score_cutoff]
        m=m[:,pa.pocketsinfo.GetAllScore()>score_cutoff]
        np.savetxt(os.path.join(out_dir,sync_name),m,fmt='%6.3f',delimiter=',')
        m=pa.pocketsinfo.PocketAcorrAnalysis()
        m=m[pa.pocketsinfo.GetAllScore()>score_cutoff]
        m=m[:,pa.pocketsinfo.GetAllScore()>score_cutoff]
        np.savetxt(os.path.join(out_dir,accor_name),m,fmt='%6.3f',delimiter=',')
        np.savetxt(os.path.join(out_dir,id_name),np.arange(0,pa.pocketsinfo.GetAllScore().shape[0])[pa.pocketsinfo.GetAllScore()>score_cutoff],delimiter=',')

    @staticmethod
    def WriteGroupVtime(pa:PocketsAnalysis,out_dir:str,**kw)->None:
        name_file=kw.get('name_file','Pocket_group_name_list.dat')
        main_file=kw.get('main_file','Pocket_group_mainname_list.dat')
        g_v_t=kw.get('v_time_file','Pocket_group_v_time.dat')
        id_name=kw.get('id_name','Pocket_sync_corr_id.dat')
        #pocket score cutoff value
        score_cutoff=kw.get('score_cutoff',10.0)
        
        name_list,pgroup_v_time=pa.pocketsinfo.GetAllPocketGroupVTime()
        with open(os.path.join(out_dir,name_file),'w') as f:
            for i in name_list:
                f.writelines(f'{i}\n')
        np.savetxt(os.path.join(out_dir,g_v_t),pgroup_v_time.transpose(),fmt='%6.3f',delimiter=',')

        size=0
        main_name=[]
        for k,v in pa.pocketsinfo._pocket_composition.items():
            size+=len(v[2])
        out=np.zeros((size,pa.size()))-1
        main_list=list(pa.pocketsinfo._pocket_composition.values())
        main_list=sorted(main_list,key=lambda x:x[1],reverse=True)
        count=0
        for i in main_list[::-1]:
            count+=len(pa.pocketsinfo._pocket_composition[i[0]][2])
            main_name.append([i[0],count])
        count=0
        with open(os.path.join(out_dir,main_file),'w') as f:
            for i in main_name:
                f.writelines(f'{i[0]},{i[1]}\n')
        np.savetxt(os.path.join(out_dir,id_name),np.arange(0,pa.pocketsinfo.GetAllScore().shape[0])[pa.pocketsinfo.GetAllScore()>score_cutoff],delimiter=',')

    @staticmethod
    def WtriteMainTransProb(pa:PocketsAnalysis,out_dir:str,**kw)->None:
        frame_cutoff=kw.get('frame_cutoff',0)
        mtp_name_file=kw.get('name_file','MTP_name_list.dat')
        mtp_matrix_file=kw.get('matrix_file','MTP_matrix.dat')
        name_list,matrix=pa.pocketsinfo.MainPocketTransProbAnalysis(frame_cutoff=frame_cutoff)
        with open(os.path.join(out_dir,mtp_name_file),'w') as f:
            for i in name_list:
                f.writelines(f'{i}\n')
        np.savetxt(os.path.join(out_dir,mtp_matrix_file),matrix)

    @staticmethod
    def WritePocketInfos(pa:PocketsAnalysis,out_dir:str,**kw)->None:
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        pock_dir=os.path.join(out_dir,'pock_data')
        if not os.path.exists(pock_dir):
            os.makedirs(pock_dir)

        np.save(os.path.join(pock_dir,'vol_v_time.npy'),pa.pocketsinfo._pocket_volume_v_time)
        np.save(os.path.join(pock_dir,'vol_v_time_sorted.npy'),pa.pocketsinfo._pocket_volume_v_time_sorted)
        np.save(os.path.join(pock_dir,'voltotal_v_time.npy'),pa.pocketsinfo.GetTotalVolVTime())
        np.save(os.path.join(pock_dir,'pock_score.npy'),pa.pocketsinfo.GetAllScore())
        np.save(os.path.join(pock_dir,'pock_lifetime.npy'),pa.pocketsinfo.GetAllLifeTime())
        np.save(os.path.join(pock_dir,'npr_v_time.npy'),pa.pocketsinfo._pocket_nonpolar_v_time)
        if pa.rec_mask.shape[0]!=0 and pa.lig_mask.shape[0]!=0:
            np.save(os.path.join(pock_dir,'occupancy_v_time.npy'),pa.pocketsinfo._pocket_occupancy_v_time)
            mean_occu=[]
            for i in range(pa.pocketsinfo.GetPocketNum()):
                mean_occu.append(pa.pocketsinfo.GetMeanOccRatio(i,80))
            np.save(os.path.join(pock_dir,'occupancy_mean.npy'),np.array(mean_occu))

def ParserBox(input:List[str])->Tuple[List[List[int]],List[List[float]]]:
    atom_id=[]
    length=[]
    tmp=[]
    for group in input:
        if len(group)!=2:
            print('Wrong number of box parameters!')
            exit()
        tmp=group[0].split(',')
        atom_id.append(list(map(int,tmp)))
        tmp=group[1].split(',')
        if len(tmp)!=3:
            print('Wrong number of box length parameters!')
            exit()
        length.append(list(map(float,tmp)))
    return (atom_id,length)

def ParserMask(protein:mdtraj.Trajectory,instr:str,mode:str='all')->List[int]:
    inlist=instr.split(',')
    inlist=list(map(int,inlist))
    resid=[]
    for i in range(0,len(inlist),2):
        resid.extend(list(range(inlist[i]-1,inlist[i+1]-1)))
    atomid=[]
    if mode=='all':
        atomid=[atom.index for atom in protein.topology.atoms if atom.residue.index in resid]
    if mode=='bb':
        atomid=[atom.index for atom in protein.topology.atoms if (atom.residue.index in resid and atom.name in ['CA','C','N'])]
    return atomid

def GetPA(md)->PocketsAnalysis:
    print('Read traj file...',end='')
    if not os.path.exists(md['top']):
            print(f'The topology file does not exist! ({md["top"]})')
            exit()
    if len(md['traj'])==0:
            print('Traj file not set!')
            exit()
    traj_list=md['traj'].split(',')
    for i in traj_list:
        i=i.strip()
        if not os.path.exists(i):
            print(f'The traj file does not exist! ({i})')
            exit()
    protein=mdtraj.load(md['traj'].split(','),top=md["top"])
    print(f'{"done":>8s}.')

    lig_mask=[]
    rec_mask=[]
    if 'lig_mask' in md.keys() and len(md['lig_mask'])!=0:
        print('Parse  lig mask...',end='')
        lig_mask=ParserMask(protein,md['lig_mask'])
        print(f'{"done":>8s}.')
    if 'rec_mask' in md.keys() and len(md['rec_mask'])!=0:
        print('Parse  rec mask...',end='')
        rec_mask=ParserMask(protein,md['rec_mask'])
        print(f'{"done":>8s}.')
    pa=PocketsAnalysis(protein,rec_mask=rec_mask,lig_mask=lig_mask)

    if 'box' in md.keys() and len(md['box'])!=0:
        print('Parse Box...',end='')
        if type(md['box'])==str:
            tmp=md['box'].split()
            if len(tmp)==0 or len(tmp)%2==1:
                print('Wrong number of box parameters!')
            box_list=[]
            for i in range(0,len(tmp),2):
                box_list.append([tmp[i],tmp[i+1]])
        else:
            box_list=md['box']
        atom_id,box_length=ParserBox(box_list)
        pa.SetBox(atom_id,box_length)
        print(f'{"done":>8s}.\n')

    return pa

def WriteFiles(pa:PocketsAnalysis,md)->None:
    out_dir=md['out']
    if bool(md.get('out_best','True')):
        print('Write best conformation...',end='')
        PAHelper.WriteBestFrames(pa,os.path.join(out_dir,'./pock_info/best_conf'),life_time_coutoff=int(md.get('best_lt_cutoff','10')),frame_shift=int(md.get('frame_start','0')),reserve_num=int(md.get('best_num','10')),is_use_scoring=bool(md.get('best_use_score','True')))
        print(f'{"done":>8s}.')
    if bool(md.get('out_pock_info','False')):
        print('Write pocket information...',end='')
        PAHelper.WritePocketInfos(pa,os.path.join(out_dir,'./pock_info/'))
        print(f'{"done":>8s}.')
    if bool(md.get('out_main_best','False')):
        print('Write main group...',end='')
        PAHelper.WriteMainBestFrames(pa,os.path.join(out_dir,'./pock_info/main_group'),life_time_coutoff=int(md.get('main_lt_cutoff','10')),frame_shift=int(md.get('frame_start','0')),reserve_num=int(md.get('main_num','1')),is_use_scoring=bool(md.get('main_use_score','True')))
        print(f'{"done":>8s}.')
    if bool(md.get('out_coex','True')):
        print('Write coexist matrix...',end='')
        PAHelper.WritePocketAccor(pa,os.path.join(out_dir,'./pock_info/'),score_cutoff=float(md.get('score_cutoff','10')))
        print(f'{"done":>8s}.')
    if bool(md.get('out_corr','True')):
        print('Write correlation matrix...',end='')
        PAHelper.WriteGroupVtime(pa,os.path.join(out_dir,'./pock_info/'))
        print(f'{"done":>8s}.')
    if bool(md.get('out_mtp','True')):
        print('Write main pocket transformation probability matrix...',end='')
        PAHelper.WtriteMainTransProb(pa,os.path.join(out_dir,'./pock_info/'),frame_cutoff=float(md.get('main_mtp_cutoff','10')))
        print(f'{"done":>8s}.\n')


#%% plot class
def cmToinch(value): 
    return value/2.54 

class PlotHalper:
    def __init__(self,**kw) -> None:
        plt.rcParams['figure.dpi'] = kw.get('figdpi',300)
        plt.rcParams['savefig.dpi'] = kw.get('sfigdpi',300)
        matplotlib.rcParams['font.family'] = kw.get('font','Times New Roman')
        self.flierprops={
        'markerfacecolor':'aliceblue',
        'marker': 'D',
        'color': 'black',
        'markeredgecolor': 'black',
        'markersize': 3.0,
        'linestyle': 'none',
        'linewidth': kw.get('linewidth',1.0)}

    def Hist(self,ay:np.ndarray,params:dict={}):
        fig, ax = plt.subplots(figsize=(cmToinch(params.get('length',20)),cmToinch(params.get('width',10))))

        ax.set_xlabel(params.get('xlabel','xlabel'))
        ax.set_ylabel(params.get('ylabel','ylabel'))
        ax.set_title(params.get('title','title'))

        ax.set_yscale(params.get('xscale','linear'))#linear log symlog logit
        ax.set_yscale(params.get('yscale','linear'))
        #ax.set_xticks(np.arange(0,40,1),[str(i) for i in range(10,51)])
        #ax.xaxis.set_major_formatter(lambda x,pos:str(pos+10))
        # if 'ylim1' in params and 'ylim2' in params:
        #     ax.set_ylim(params['ylim1'],params['ylim2'])

        ax.hist(ay,bins=params.get('bins',10),range=params.get('range',None),density=params.get('density',False),cumulative=params.get('cumulative',False))
        #plt.xticks(rotation=params.get('xticks_rotation',0))
        #plt.yticks(rotation=params.get('yticks_rotation',0))
        return fig, ax

    def test(self,ay:np.ndarray,params:dict={},axes=None,ax_row=-1,ax_col=-1):
        ax=axes
        ax[ax_row,ax_col].set_xlabel(params.get('xlabel','xlabel'))
        ax[ax_row,ax_col].set_ylabel(params.get('ylabel','ylabel'))
        ax[ax_row,ax_col].set_title(params.get('title','title'))

        ax[ax_row,ax_col].set_yscale(params.get('xscale','linear'))#linear log symlog logit
        ax[ax_row,ax_col].set_yscale(params.get('yscale','linear'))

        ax[ax_row,ax_col].hist(ay,bins=params.get('bins',10),range=params.get('range',None),density=params.get('density',False),cumulative=params.get('cumulative',False))
#%%arg parser
class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text: str, width: int):
        if text.startswith('R|'):
            return text[2:].splitlines()
        return super()._split_lines(text, width)

parser = argparse.ArgumentParser(
                    prog='PocketAnalysis',
                    formatter_class=SmartFormatter,
                    description='Analyzing protein pocket features in molecular dynamics simulations trajectories.')

parser.add_argument('--version', action='version', version='%(prog)s 1.1.2')
parser.add_argument('--top', type=str,metavar='Path_to_your_topology_file',default='',
                    help='This parameter specifies the path to your topology file.')
parser.add_argument('--traj', type=str,metavar='Path_to_your_traj_file',default='',
                    help="This parameter specifies the path to the tracking file. Multiple paths separated by ','")
parser.add_argument('--rec_mask',type=int,action='extend',nargs='+',default=[],metavar='rec_id',
                    help='''R|This parameter specifies your receptor protein residue index. 
If neither rec_mask nor lig_mask are specified, all atoms are considered as proteins. 
If rec_mask is not specified but lig_mask is specified, then rec_mask = total_atom - lig_mask.
(Note: Parameters must appear in pairs.) 
Format: --rec_mask resn1 resn2 resn3 resn4... 
Example: If you want to specify residues 1-50, 60 and 70-89 as your receptor protein: "-rec_mask 1 51 60 61 70 90"
''')
parser.add_argument('--lig_mask',type=str,default='',metavar='lig_atom_id',
                    help='This parameter specifies the atomic index of your ligand. The default is null. See "--rec_mask" for parameter format')

parser.add_argument('--box',type=str,action='append',nargs='+',default=[],metavar='x',
                    help='''R|This parameter specifies the position boxes.
If you are only interested in a certain region of the protein, you can specify one or more location boxes that wrap around that region. 
The analysis will be focused within the location boxes, which allows the program output to filter out data you are not interested in.
At the same time, this can speed up the running of the program and reduce the memory requirements.
NOTE: To accommodate continuous changes in protein conformation, the box centers are set to the coordinates of a single atom or to the geometric centers of several atoms, 
which allows the box coordinates to be synchronized with the conformational changes of the protein.
Format: --box atom1_id,atom2_id length(X-axis direction),width(Y-axis direction),height(Z-axis direction)
Example: --box 372 18,10,22 --box 458,963 14,12,20''')

parser.add_argument('--dist_cutoff',type=float,default=3.0,metavar='float, default=3.0',help='Distance cutoff value when clustering subpockets')


parser.add_argument('--frame_start',type=int,default=0,metavar='integer, default=0',help='This parameter specifies which frame of the trajectory to start processing from, and frames before this parameter will be ignored')
parser.add_argument('--frame_stop',type=int,default=-1,metavar='integer, default=-1',help='This parameter specifies the frame to which the program ends processing, and the frames after this parameter will be ignored. The default value of this parameter is the last frame')
parser.add_argument('--frame_offset',type=int,default=1,metavar='integer, default=1',help='This parameter specifies the frame interval step size for processing trajectories.')

parser.add_argument('--out',type=str,metavar='out_dir',default='./', help='The folder where the analysis results are saved')
parser.add_argument('--out_best',type=bool,metavar='bool, default=True',default=True,help='Output the optimal conformation, which means the overall score of the pocket is the highest')
parser.add_argument('--out_pock_info',type=bool,metavar='bool, default=False',default=False,help='Save the original information of the pocket in. npy format. For further analysis')
parser.add_argument('--out_main_best',type=bool,metavar='bool, default=False',default=False,help='Output the optimal conformation of each main_group')
parser.add_argument('--out_coex',type=bool,metavar='bool, default=False',default=False,help='Control whether to output sub pocket coexistence matrix')
parser.add_argument('--out_corr',type=bool,metavar='bool, default=False',default=False,help='Control whether to output sub pocket correlation matrix')
parser.add_argument('--out_mtp',type=bool,metavar='bool, default=False',default=False,help='Control whether to output main pocket transformation probability matrix')

parser.add_argument('--best_num',type=int,default=10,metavar='integer, default=10',help='Maximum number of best conformation retention')
parser.add_argument('--best_lt_cutoff',type=int,default=10,help='Truncated value of lifetime. Pockets with a lifetime less than this time will be ignored when calculating pocket scoring.')
parser.add_argument('--best_use_score',type=bool,default=True,metavar='integer, default=10, unit: frame',help='Whether to use pocket scoring function when calculating the optimal conformation. If it is not, the pocket volume will be used to rate the pocket.')
parser.add_argument('--main_num',type=int,default=1,help='Maximum number of best main group conformation retention.')
parser.add_argument('--main_lt_cutoff',type=int,default=10,help='Truncated value of lifetime. Pockets with a lifetime less than this time will be ignored when calculating pocket scoring.')
parser.add_argument('--main_use_score',type=bool,default=True,metavar='bool, default=True',help='Whether to use pocket scoring function when calculating the optimal main group conformation.')
parser.add_argument('--main_mtp_cutoff',type=int,default=0,metavar='integer, default=0',help='Truncated value of main group lifetime.')

parser.add_argument('--config',type=str,default='',help='Specify the path to the control file. When this parameter is specified, all input information will be read from the control file. Other command line input parameters will be ignored')
#%%
# a=' --traj ./traj.nc --lig_mask 612 613 --frame_start 500 --frame_stop 1000 --out ./test --out_pock_info true --out_main_best true --out_coex true --out_corr true'
# md=vars(parser.parse_args(a.split()))
#(pytorch) PS G:\data\sars-cov-2\main_protease\mp_ihb_ensitrelvir\analysis\alphaspace> python G:\data\pocket_analysis\project\pocket_analysis.py --top .\com_wat_strip.prmtop --traj .\traj.nc --lig_mask 612 613 --frame_start 500 --frame_stop 1000 --out ./test --out_pock_info true --out_main_best true --out_coex true --out_corr true
#conda init bash
#%%
if __name__=='__main__':
    if len(sys.argv)==1:
        parser.parse_args(['-h'])
    args=parser.parse_args(sys.argv[1:])
    config={}
    
    config_path=args.config
    if args.config!='':
        print('Control file specified, attempting to load parameters from control file.')
        if not os.path.exists(config_path):
            print('The config file does not exist!')
            exit()
        config = ConfigParser()
        config.read(args.config)
    else:
        config={'GENERAL':vars(args),'MODEL0':vars(args)}
        config['GENERAL']['mode']='single'

    if config['GENERAL']['mode']=='single':
        print('single mode')
        pa=GetPA(config['MODEL0'])
        pa.SetDistCutoff(float(config['GENERAL'].get('dist_cutoff','3.0')))
        ts=int(config[f'MODEL0'].get('frame_start','0'))
        te=int(config[f'MODEL0'].get('frame_stop','-1'))
        to=int(config[f'MODEL0'].get('frame_offset','1'))
        pa.Analysis(ts,te,to)
        WriteFiles(pa,config['MODEL0'])
        for i in range(10):
            pa.rec[i].save(f'./test{i}.pdb')
    elif config['GENERAL']['mode']=='multi':
        print('multi mode')
        models=ModelGroup()
        models.SetDistCutoff(float(config['GENERAL'].get('dist_cutoff','3.0')))
        align_id=[]
        step=[]

        for i in range(100):
            if f'MODEL{i}' in config:
                print(f'MODEL{i}')
                # for k,v in config[f'MODEL{i}'].items():
                #     print(f'{k} : {v}')
                models.AddPA(GetPA(config[f'MODEL{i}']))
                align_id.append(ParserMask(models[i].rec,config[f'MODEL{i}'].get('align_mask').strip(),mode='bb'))
                ts=int(config[f'MODEL{i}'].get('frame_start','0'))
                te=int(config[f'MODEL{i}'].get('frame_stop','-1'))
                to=int(config[f'MODEL{i}'].get('frame_offset','1'))
                step.append([ts,te,to])
        print('start align model...    ',end='')
        models.AlignModel(align_id)  
        print('done')
        models.Analysis(step)

        for i in range(len(models.pa_list)):
            print(f'Write data of MODEL{i}')
            WriteFiles(models.pa_list[i],config[f'MODEL{i}'])
