#%%
import os
import sys
import argparse
import pickle
import alphaspace2 as al
from alphaspace2 import Snapshot as ss
from typing import List,Union,Dict,Tuple
import mdtraj
import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.linalg import orthogonal_procrustes
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
        
        self.lig_cutoff=3.0

    def __getitem__(self,id:int):
        return self.sslist[id]

    def SetBox(self,center_atom_id:List[List[int]],length:List[List[float]]):
        self.boxCenterAtomId=center_atom_id
        self.boxlen=length

    def SetLigCutoff(self,lig_cutoff:float):
        self.lig_cutoff=lig_cutoff

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
            self.sslist.append(al.Snapshot(lig_cutoff=self.lig_cutoff))
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
        self._pocket_rank:np.ndarray
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
        #self._pocket_xyz:Dict[int,np.ndarray]={}#ingore invalid sub-pocket,{pid,n*3 array}
        self._total_score_v_time:np.ndarray
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
                main_key.append([self._pocket_rank[j]+1,self.GetScore(j)])
                main_score+=self.GetScore(j)
            elif self.GetScore(j)>sub_score_cutoff:
                sub_key.append([self._pocket_rank[j]+1,self.GetScore(j)])
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

    def _CalcTotalScoreVTime(self,score_cutoff:float=0,percentile:int=80,is_use_scoring:bool=True):
        self._total_score_v_time=np.zeros(self.GetFrameNum())
        score=self._pockets_score.copy()
        score[score<score_cutoff]=0
        score[score>0]=np.log10(score[score>0]-8)
        if is_use_scoring==False:
            score=score*0+1
        for i in range(self.GetFrameNum()):
            tsumscore=0
            for j in range(self.GetPocketNum()):
                    tvol=self._pocket_volume_v_time[j,i] if self._pocket_volume_v_time[j,i]<self.GetPCTLVol(j,percentile) else self.GetPCTLVol(j,percentile)
                    tsumscore+=(tvol*score[j])
            self._total_score_v_time[i]=tsumscore

    def _GetBestConf(self,frame_list:List[int],reserve_num:int)->list:
        ttscore=self._total_score_v_time.copy()[frame_list]
        tvol=self.GetTotalVolVTime(score_cutoff=10)
        sort_id=np.argsort(ttscore)
        reserve_num=reserve_num if reserve_num<len(frame_list) else len(frame_list)
        out=[]
        for i in range(reserve_num):
            tid=-1-i
            fid=frame_list[sort_id[tid]]
            out.append([fid,self._total_score_v_time[fid],tvol[fid],self.GetMainName(fid),self.GetMinorName(i)])
        return out

    def _CalaPocketRank(self)->None:
        all_p_rank=[[i,self.GetScore(i)] for i in range(self.GetPocketNum()) ]
        all_p_rank=np.array(sorted(all_p_rank,key=lambda x:x[1],reverse=True))[:,0].transpose()
        self._pocket_rank=np.zeros(self._pockets_score.shape[0],dtype=int)
        for i in range(len(all_p_rank)):
            self._pocket_rank[int(all_p_rank[i])]=i

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
        '''
        percentage:int 0-100
        '''
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

    def GetTotalVolVTime(self,score_cutoff=0)->np.ndarray:
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
        if l>10:
            tnpr=tnpr[int((100-percentile)/200*l):int(np.round((0.5+percentile/200)*l))]
        return np.mean(np.power(tnpr-np.mean(tnpr),2))

    def GetNPRVTime(self,pid:int)->np.ndarray:#0==nan
        return self._pocket_nonpolar_v_time[pid].copy()

    def GetNoZeroNPRVTime(self,pid:int)->np.ndarray:
        tnpr=self._pocket_nonpolar_v_time[pid]
        tnpr=tnpr[tnpr>0]
        return tnpr.copy()

    def GetNPRFluctuations(self,pid:int,percentage:int)->Tuple[float,float,float,float]:
        tnpr=self._pocket_nonpolar_v_time[pid]
        tnpr=np.sort(tnpr)
        tnpr=tnpr[self._pocket_volume_v_time_sorted[pid]>0]
        if tnpr.shape[0]==0:
            return (0,0,0,0)
        if tnpr.shape[0]==1:
            return (tnpr[0],tnpr[0],0,0)
        if tnpr.shape[0]<10:
            tmin=tnpr[0]
            tmax=tnpr[-1]
        else:
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
    
    def GetRank(self,pid:int)->float:
        return self._pocket_rank[pid]

    def GetAllScore(self)->np.ndarray:
        return self._pockets_score.copy()

    def GetAllRank(self)->np.ndarray:
        return self._pocket_rank.copy()

    def PocketCoexAnalysis(self,score_cutoff:float=.0)->Tuple[np.ndarray,np.ndarray,np.ndarray]:
        pnum=self.GetPocketNum()
        out=np.eye(pnum)
        pid=np.arange(pnum,dtype=int)
        prank=self._pocket_rank.copy().astype(int)
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
        if score_cutoff>0:
            out=out[self.GetAllScore()>score_cutoff]
            out=out[:,self.GetAllScore()>score_cutoff]
            pid=pid[self.GetAllScore()>score_cutoff]
            prank=prank[self.GetAllScore()>score_cutoff]
        return (pid,prank,out)

    def PocketAcorrAnalysis(self,score_cutoff:float=.0)->Tuple[np.ndarray,np.ndarray,np.ndarray]:
        pnum=self.GetPocketNum()
        out=np.eye(pnum)
        pid=np.arange(pnum,dtype=int)
        prank=self._pocket_rank.copy().astype(int)
        for i in range(pnum):
            for j in range(i+1,pnum):
                pivt=(self.GetVolVTime(i)[1::]-self.GetVolVTime(i)[:-1:])/(self.GetMaxVol(i)+1e-5)
                povt=(self.GetVolVTime(j)[1::]-self.GetVolVTime(j)[:-1:])/(self.GetMaxVol(j)+1e-5)
                out[i,j]=out[j,i]=np.sum(pivt*povt)/(np.linalg.norm(pivt)*np.linalg.norm(povt)+1e-5)
        if score_cutoff>0:
            out=out[self.GetAllScore()>score_cutoff]
            out=out[:,self.GetAllScore()>score_cutoff]
            pid=pid[self.GetAllScore()>score_cutoff]
            prank=prank[self.GetAllScore()>score_cutoff]
        return (pid,prank,out)

    def GetBestConf(self,reserve_num:int,percentile:int=80,is_use_scoring:bool=True):
        frame_list=list(range(self.GetFrameNum()))
        return self._GetBestConf(frame_list,reserve_num)

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

    def GetMainBestConf(self,main_name:str,reserve_num:int):
        return self._GetBestConf(self._pocket_composition[main_name][3],reserve_num)

    def GetGroupBestConf(self,main_name:str,sub_name:str,reserve_num:int):
        return self._GetBestConf(self._pocket_composition[main_name][2][sub_name][2],reserve_num,)

    def MainPocketTransProbAnalysis(self)->Tuple[List[str],List[int],np.ndarray,np.ndarray]:
        main_list=list(self._pocket_composition.values())
        main_list=sorted(main_list,key=lambda x:x[1],reverse=True)#sort by score
        main_frames_num_list=np.array([len(i[-1]) for i in main_list])
        mname_list=[i[0] for i in main_list]
        out=np.zeros((len(main_list),len(main_list)))
        m_v_t=self.GetMainGroupVTime()
        for i in range(len(m_v_t)-1):
            #if m_v_t[i]!=m_v_t[i+1]:
            out[mname_list.index(m_v_t[i]),mname_list.index(m_v_t[i+1])]+=1
        mname_list=[i[0] for i in main_list]
        mlt_list=[len(i[-1]) for i in main_list]
        norm_out=out.copy()
        for i in range(len(main_list)):
            norm_out[i]/=len(main_list[i][-1])
        return (mname_list,mlt_list,out,norm_out)

class PocketsAnalysis:
    def __init__(self,rec:mdtraj.Trajectory,rec_mask:List[int]=[],lig_mask:List[int]=[],**kw) -> None:
        self.rec:mdtraj.Trajectory=rec
        self.lig_mask:np.ndarray=np.array(lig_mask,dtype=int)
        self.rec_mask:np.ndarray=np.array(rec_mask,dtype=int)
        if self.lig_mask.shape[0]!=0:
            self.lig_mask-=1
        if self.rec_mask.shape[0]==0:
            self.rec_mask:np.ndarray=np.arange(0,self.rec.xyz.shape[1])
            self.rec_mask=np.setdiff1d(self.rec_mask,self.lig_mask,True)
        else:
            self.rec_mask-=1
        
        self.rec.superpose(self.rec,0,self.rec_mask)
        self.snap_shots=SnapShots(rec)
        self.pocketsinfo=PocketsInfo()
        self.distance_cutoff:float=3.0

        self.score_cutoff=kw.get('score_cutoff',10.0)
        self.is_use_score=kw.get('is_use_score',True)
        self.percentile=kw.get('percentile',80.0)
        #score_cutoff:float=0,percentile:int=80,is_use_scoring:bool=True

    def size(self)->int:
        '''
        return frame number
        '''
        return len(self.snap_shots.sslist)
    
    def _ClusterPockets(self):
        total_coords=[]
        #Extract pocket coordinates
        for ss in self.snap_shots.sslist:
            total_coords.append(ss._pocket_xyz)
        total_coords=np.concatenate(total_coords,axis=0)
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
        
    def GetPocketXYZ(self,pid:int)->np.ndarray:
        out=[]
        for i in range(len(self.snap_shots.sslist)):
            tpids=self.pocketsinfo._pocket_cluster_splited[i]#tmp pockets index tpids:pock index list
            if pid in tpids:
                tpids_count=np.bincount(tpids)
                parg=np.argwhere(tpids==pid).flatten()
                tcoord=np.zeros((1,3),dtype=float)
                for j in parg:
                    tcoord+=self.snap_shots[i]._pocket_xyz[j]
                tcoord/=tpids_count[pid]
                out.append(tcoord)
        return np.concatenate(out,axis=0)

    def _rotateMatrix(self,M:np.ndarray,N:np.ndarray)->np.ndarray:
        center_M = M - np.mean(M, axis=0)
        center_N = N - np.mean(N, axis=0)
        R, _ = orthogonal_procrustes(center_M, center_N)
        return R

    def Align(self,coord:np.ndarray,mask:np.ndarray):
        if coord.shape[0]!=mask.shape[0]:
            print('Atom number mismatch! Unable to perform pocket matching, skip.',end='')
            return False
        N=self.rec.xyz[0,mask,:]
        R=self._rotateMatrix(N,coord)
        T=-np.mean(N,axis=0)
        Ta=T*10
        for i in range(self.size()):
            self.rec.xyz[i]+=T
            self.rec.xyz[i]=self.rec.xyz[i]@R
            self.snap_shots[i]._alpha_xyz+=Ta
            self.snap_shots[i]._alpha_xyz=self.snap_shots[i]._alpha_xyz@R
            self.snap_shots[i]._pocket_xyz+=Ta
            self.snap_shots[i]._pocket_xyz=self.snap_shots[i]._pocket_xyz@R
            #self.snap_shots[i].Translation(T)
            #self.snap_shots[i].Rotate(R)
        return True

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

    def isUseScore(self,is_use_score:bool):
        self.is_use_score=is_use_score

    def SetLigCutoff(self,lig_cutoff:float):
        self.snap_shots.SetLigCutoff(lig_cutoff)

    def SetDistCutoff(self,distance:float):
        self.distance_cutoff=distance

    def SetScoreCutoff(self,score_cutoff:float):
        self.score_cutoff=score_cutoff

    def SetPercentile(self,percentile:int):
        if percentile<0 or percentile>100:
            print('Percentile must be an integer between 0 and 100')
            return
        self.percentile=percentile

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
        self.pocketsinfo._CalaPocketRank()
        self.pocketsinfo._CalcTotalScoreVTime(score_cutoff=self.score_cutoff,percentile=self.percentile,is_use_scoring=self.is_use_score)
        self.pocketsinfo._PocketGrouppAnalysis()
        print(' done.\n')

    def WritePockets(self,frame:int,out_file:str)->None:
        resi=1
        serial=1
        piday=self.pocketsinfo._pocket_cluster_splited[frame]
        prank=[self.pocketsinfo._pocket_rank[i] for i in piday]
        #sort_id=np.argsort(prank)
        pxyz=self.snap_shots.sslist[frame]._pocket_xyz
        axyz=self.snap_shots.sslist[frame]._alpha_xyz
        #print(f'piday={piday}')
        with open(out_file,'w') as f:
            tline=''
            for indx in range(piday.shape[0]):
                #print(f'indx={indx},pockcomp={self.snap_shots.GetPockComposition(frame,indx)}')
                #serial[]atom_name[]res_name[][2]res_seq[][][][4]xyz-occupancy-tempfactor[][][][][][][][][][][11]element-charge
                tline=f"{'ATOM':<6s}{serial:>5d} {'PCC':^4s} {'APL':>3s}  {resi:>4d}    {pxyz[indx][0]:8.3f}{pxyz[indx][1]:8.3f}{pxyz[indx][2]:8.3f}{1.0:6.2f}{self.pocketsinfo._pocket_rank[piday[indx]]+1:6.2f}          {'PC':2s} 0\n"
                f.writelines(tline)
                serial+=1
                for aid in self.snap_shots.GetPockComposition(frame,indx):
                    #print(f'aid={aid}')
                    if self.snap_shots.sslist[frame]._alpha_contact[aid]==0:
                        tline=f"{'ATOM':<6s}{serial:>5d} {'AAU':^4s} {'APL':>3s}  {resi:>4d}    {axyz[aid,0]:8.3f}{axyz[aid,1]:8.3f}{axyz[aid,2]:8.3f}{1.0:6.2f}{self.pocketsinfo._pocket_rank[piday[indx]]+1:6.2f}          {'AC':2s} 0\n"
                    else:
                        tline=f"{'ATOM':<6s}{serial:>5d} {'AAO':^4s} {'APL':>3s}  {resi:>4d}    {axyz[aid,0]:8.3f}{axyz[aid,1]:8.3f}{axyz[aid,2]:8.3f}{1.0:6.2f}{self.pocketsinfo._pocket_rank[piday[indx]]+1:6.2f}          {'AC':2s} 0\n"
                    f.writelines(tline)
                    serial+=1
                resi+=1

class ModelGroup:
    def __init__(self,**kw) -> None:
        self.pa_list:List[PocketsAnalysis]=[]
        self.need_analysis:List[int]=[]
        self.distance_cutoff=kw.get('disCutoff',3.0)

    def __getitem__(self,id:int)->PocketsAnalysis:
        return self.pa_list[id]

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

    def AlignModel(self,align_mask:List[np.ndarray])->None:
        '''
        align_mask:[ [system1 atom indexs],   [system2 atom indexs] .....   ]
        '''
        #center model 1
        coord_ref=self.pa_list[0].rec.xyz[0,align_mask[0],:]
        T=np.mean(coord_ref,axis=0)
        Ta=T*10
        coord_ref=self.pa_list[0].rec.xyz[0,align_mask[0],:]
        for i in range(self.pa_list[0].size()):
            self.pa_list[0].rec.xyz[i]-=T
            self.pa_list[0].snap_shots[i]._alpha_xyz-=Ta
            self.pa_list[0].snap_shots[i]._pocket_xyz-=Ta
        #align model 2...
        for i in range(1,self.size()):
            print(f'Align model{i}...',end='')
            self.pa_list[i].Align(coord_ref,align_mask[i])
            print(f'{"done":>8s}.')

    def PocketMatch(self,score_cutoff:float=0.0,dist_cutoff=3.0):
        tmcoord=[]
        tmpid=[]
        tmsplit=[]
        for pa in self.pa_list:
            pa_pscore=pa.pocketsinfo.GetAllScore()
            pa_pid=np.argwhere(pa_pscore>score_cutoff).flatten()
            tmsplit.append(pa_pid.shape[0])
            for pid in pa_pid:
                pcoord=np.mean(pa.GetPocketXYZ(pid),axis=0)
                tmcoord.append(pcoord)
                tmpid.append(pid)
        tmcoord=np.stack(tmcoord,axis=0)
        cid=fcluster(linkage(tmcoord, method='average'), dist_cutoff,criterion='distance')
        cid-=1
        return (np.array(tmsplit),np.array(tmpid),cid)

    def Analysis(self,step:List[List[int]])->None:
        for i in self.need_analysis:
            print(f'Analysis MODEL{i}....     ')
            self.pa_list[i].Analysis(step[i][0],step[i][1],step[i][2])
            print('Analysis done.\n')

#%% helper class
class PAHelper:
    @staticmethod
    def _WriteBestFrames(pa:PocketsAnalysis,out_dir:str,result:list,**kw)->None:
        frame_shift=kw.get('frame_shift',0)
        is_out_vol_v_time=kw.get('vol_v_time',True)
    
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

        pock_rank=pa.pocketsinfo.GetAllRank()

        out=result
        info=[]
        for i in range(len(out)):
            pcom=pa.pocketsinfo.GetPocketComposition(out[i][0])
            tmp=[[pcom[j],pa.pocketsinfo.GetScore(pcom[j]),pa.pocketsinfo.GetVol(pcom[j],out[i][0])]for j in range(len(pcom))]
            tmp=sorted(tmp,key=lambda x:x[1],reverse=True)
            info.append(tmp)
            pa.WritePockets(out[i][0],os.path.join(pdb_dir,f'apocket-r{i+1:03d}.pdb'))
            pa.rec[frame_shift+out[i][0]].save(os.path.join(pdb_dir,f'protein-r{i+1:03d}-f{frame_shift+out[i][0]*pa.offset:05d}.pdb'))
        with open(os.path.join(out_dir,'out.dat'),'w') as f:
            for i in range(len(out)):
                if out[i][0]==-1:
                    continue
                f.writelines(f'rank={i+1} frame={out[i][0]*pa.offset+frame_shift:6d} pock_score={out[i][1]:10.3f} pock_space={out[i][2]} main_name={out[i][3]} minor_name={out[i][4]}\n{"id":^3s} {"rank":^4s} {"socre":^6s} {"space":^7s} {"mnpr%80":^8s} {"occu":^6s} {"m_occu":^6s}\n')
                for j in info[i]:
                    f.writelines(f'{j[0]:3d} {pock_rank[j[0]]+1:4d} {j[1]:6.2f} {j[2]:7.2f} {pa.pocketsinfo.GetMeanNPR(j[0],80)*100:8.2f} {pa.pocketsinfo.GetOccupancyRatio(j[0],out[i][0]):6.4f} {mean_occu[j[0]]:6.4f}\n')
                f.writelines('\n\n')

        if is_out_vol_v_time:
            np.savetxt(os.path.join(out_dir,'vol_effect_v_time.dat'),pa.pocketsinfo.GetTotalVolVTime(score_cutoff=10))
            np.savetxt(os.path.join(out_dir,'vol_total_v_time.dat'),pa.pocketsinfo.GetTotalVolVTime(score_cutoff=0))
        
    @staticmethod
    def WriteBestFrames(pa:PocketsAnalysis,out_dir:str,**kw)->None:
        '''
        reserve_num(type:int default:10)
        frame_shift(type:int default:0)
        '''
        reserve_num=kw.get('reserve_num',10)
        frame_shift=kw.get('frame_shift',0)
        #score_cutoff=kw.get('score_cutoff',10)

        out=pa.pocketsinfo.GetBestConf(reserve_num=reserve_num)
        PAHelper._WriteBestFrames(pa,out_dir,out,frame_shift=frame_shift)

    @staticmethod
    def WriteMainBestFrames(pa:PocketsAnalysis,out_dir:str,**kw)->None:
        '''
        reserve_num(type:int default:10)
        frame_shift(type:int default:0)
        '''
        reserve_num=kw.get('reserve_num',1)
        percentile=kw.get('percentile',80)
        is_use_scoring=kw.get('is_use_scoring',True)
        frame_shift=kw.get('frame_shift',0)

        main_list=list(pa.pocketsinfo._pocket_composition.values())
        main_list=sorted(main_list,key=lambda x:x[1],reverse=True)
        for i in main_list:
            m_out_dir=os.path.join(out_dir,f'main_{i[0]}_s{np.round(i[1])}')
            data=pa.pocketsinfo.GetMainBestConf(i[0],reserve_num=reserve_num)
            PAHelper._WriteBestFrames(pa,m_out_dir,data,frame_shift=frame_shift,vol_v_time=False)

    @staticmethod
    def WritePocketAccor(pa:PocketsAnalysis,out_dir:str,**kw)->None:
        score_cutoff=kw.get('score_cutoff',10.0)
        coex_sort_by=kw.get('coex_sort_by','id')
        corr_sort_by=kw.get('corr_sort_by','id')
        sync_name=kw.get('coexist_name','Pocket_coexistence rate.dat')
        accor_name=kw.get('accor_name','Pocket_correlation.dat')
        coex_id_name=kw.get('coex_id_name','Pocket_coex_id.dat')
        corr_id_name=kw.get('corr_id_name','Pocket_corr_id.dat')

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        pid1,prank1,m1=pa.pocketsinfo.PocketCoexAnalysis(score_cutoff=score_cutoff)
        pid2,prank2,m2=pa.pocketsinfo.PocketAcorrAnalysis(score_cutoff=score_cutoff)

        if coex_sort_by=='rank':
            sort1=np.argsort(prank1)
            pid1=pid1[sort1]
            prank1=prank1[sort1]
            m1=m1[sort1]
            m1=m1[:,sort1]
        if corr_sort_by=='rank':
            sort2=np.argsort(prank2)
            pid2=pid2[sort2]
            prank2=prank2[sort2]
            m2=m2[sort2]
            m2=m2[:,sort2]

        np.savetxt(os.path.join(out_dir,sync_name),m1,fmt='%6.3f',delimiter=',')
        np.savetxt(os.path.join(out_dir,accor_name),m2,fmt='%6.3f',delimiter=',')
        np.savetxt(os.path.join(out_dir,coex_id_name),np.stack([pid1,prank1],axis=-1),delimiter=',',header='id,rank')
        np.savetxt(os.path.join(out_dir,corr_id_name),np.stack([pid2,prank2],axis=-1),delimiter=',',header='id,rank')

    @staticmethod
    def WriteGroupVtime(pa:PocketsAnalysis,out_dir:str,**kw)->None:
        name_file=kw.get('name_file','Pocket_group_name_list.dat')
        main_file=kw.get('main_file','Pocket_group_mainname_list.dat')
        g_v_t=kw.get('v_time_file','Pocket_group_v_time.dat')
        mg_v_t_f=kw.get('mg_v_time_file','Pocket_maingroup_v_time.dat')
        #pocket score cutoff value
        score_cutoff=kw.get('score_cutoff',10.0)
        
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

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
        main_list=sorted(main_list,key=lambda x:x[1])
        mg_v_t=np.zeros((len(main_list),pa.size()))
        for i in range(len(main_list)):
            mg_v_t[i][main_list[i][-1]]=i+1
        mg_v_t-=1
        np.savetxt(os.path.join(out_dir,mg_v_t_f),mg_v_t.transpose(),fmt='%6.3f',delimiter=',')
        for i in main_list[::-1]:
            main_name.append([i[0],len(pa.pocketsinfo._pocket_composition[i[0]][-1])])
        with open(os.path.join(out_dir,main_file),'w') as f:
            f.writelines('#main_name,lifetime\n')
            for i in main_name:
                f.writelines(f'{i[0]},{i[1]}\n')

    @staticmethod
    def WtriteMainTransProb(pa:PocketsAnalysis,out_dir:str,**kw)->None:
        mtp_name_file=kw.get('name_file','MTP_name_list.dat')
        mtp_matrix_file=kw.get('matrix_file','MTP_matrix.dat')
        mtp_norm_matrix_file=kw.get('matrix_file','MTP_norm_matrix.dat')

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        name_list,life_time,matrix,norm_m=pa.pocketsinfo.MainPocketTransProbAnalysis()
        with open(os.path.join(out_dir,mtp_name_file),'w') as f:
            for i in range(len(name_list)):
                f.writelines(f'{name_list[i]} {life_time[i]}\n')
        np.savetxt(os.path.join(out_dir,mtp_matrix_file),matrix)
        np.savetxt(os.path.join(out_dir,mtp_norm_matrix_file),norm_m)

    @staticmethod
    def WritePocketInfos(pa:PocketsAnalysis,out_dir:str,**kw)->None:
        ftype=kw.get('ftype','npy')
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        pock_dir=os.path.join(out_dir,'pock_data')
        if not os.path.exists(pock_dir):
            os.makedirs(pock_dir)
        if ftype=='npy':
            np.save(os.path.join(pock_dir,'vol_v_time.npy'),pa.pocketsinfo._pocket_volume_v_time)
            np.save(os.path.join(pock_dir,'vol_v_time_sorted.npy'),pa.pocketsinfo._pocket_volume_v_time_sorted)
            np.save(os.path.join(pock_dir,'voltotal_v_time.npy'),pa.pocketsinfo.GetTotalVolVTime())
            np.save(os.path.join(pock_dir,'pock_score.npy'),pa.pocketsinfo.GetAllScore())
            np.save(os.path.join(pock_dir,'pock_lifetime.npy'),pa.pocketsinfo.GetAllLifeTime())
            np.save(os.path.join(pock_dir,'npr_v_time.npy'),pa.pocketsinfo._pocket_nonpolar_v_time)
            np.save(os.path.join(pock_dir,'scoretotal_v_time.npy'),pa.pocketsinfo._total_score_v_time)
            if pa.rec_mask.shape[0]!=0 and pa.lig_mask.shape[0]!=0:
                np.save(os.path.join(pock_dir,'occupancy_v_time.npy'),pa.pocketsinfo._pocket_occupancy_v_time)
                mean_occu=[]
                for i in range(pa.pocketsinfo.GetPocketNum()):
                    mean_occu.append(pa.pocketsinfo.GetMeanOccRatio(i,80))
                np.save(os.path.join(pock_dir,'occupancy_mean.npy'),np.array(mean_occu))
        elif ftype=='txt':
            np.savetxt(os.path.join(pock_dir,'vol_v_time.txt'),pa.pocketsinfo._pocket_volume_v_time)
            np.savetxt(os.path.join(pock_dir,'vol_v_time_sorted.txt'),pa.pocketsinfo._pocket_volume_v_time_sorted)
            np.savetxt(os.path.join(pock_dir,'voltotal_v_time.txt'),pa.pocketsinfo.GetTotalVolVTime())
            np.savetxt(os.path.join(pock_dir,'pock_score.txt'),pa.pocketsinfo.GetAllScore())
            np.savetxt(os.path.join(pock_dir,'pock_lifetime.txt'),pa.pocketsinfo.GetAllLifeTime())
            np.savetxt(os.path.join(pock_dir,'npr_v_time.txt'),pa.pocketsinfo._pocket_nonpolar_v_time)
            np.savetxt(os.path.join(pock_dir,'scoretotal_v_time.txt'),pa.pocketsinfo._total_score_v_time)
            if pa.rec_mask.shape[0]!=0 and pa.lig_mask.shape[0]!=0:
                np.savetxt(os.path.join(pock_dir,'occupancy_v_time.txt'),pa.pocketsinfo._pocket_occupancy_v_time)
                mean_occu=[]
                for i in range(pa.pocketsinfo.GetPocketNum()):
                    mean_occu.append(pa.pocketsinfo.GetMeanOccRatio(i,80))
                np.savetxt(os.path.join(pock_dir,'occupancy_mean.txt'),np.array(mean_occu))

    @staticmethod
    def WriteSummaryInfo(pa:PocketsAnalysis,out_dir:str,**kw)->None:
        fname=kw.get('fname','summary.dat')

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        with open(os.path.join(out_dir,fname),'w') as f:
            pscore=pa.pocketsinfo.GetAllScore()
            prank=pa.pocketsinfo.GetAllRank()
            pmain_num=pscore[pscore>=20].shape[0]
            pminor_num=pscore[pscore>=10].shape[0]-pscore[pscore>=20].shape[0]
            pinvalid_num=pscore[pscore<10].shape[0]
            f.writelines(f"sub-pocket num: {pa.pocketsinfo.GetPocketNum()}\nmain sub-pocket num: {pmain_num}\nminor sub-pocket num: {pminor_num}\ninvalid sub-pocket num: {pinvalid_num}\n\n")
            f.writelines(f"{'id':>5s} {'rank':>5s} {'score':>10s} {'lifetime':>10s} {'med_vol':>10s} {'max_vol':>10s} {'min_vol':>10s} {'vol_f0.8':>10s} {'moccu':>10s} {'mnpr':>10s} {'npr-std':>10s}\n")
            for i in range(pa.pocketsinfo.GetPocketNum()):
                f.writelines(f"{i:>5d} {prank[i]+1:>5d} {pscore[i]:>10.4f} {pa.pocketsinfo.GetLifeTime(i):>10d} {pa.pocketsinfo.GetMedianVol(i):>10.4f} {pa.pocketsinfo.GetMaxVol(i):>10.4f} {pa.pocketsinfo.GetMinVol(i):>10.4f} {pa.pocketsinfo.GetVolumeFluctuations(i,80)[-1]:>10.4F} {pa.pocketsinfo.GetMeanOccRatio(i):>10.4f} {pa.pocketsinfo.GetMeanNPR(i):>10.4f} {np.sqrt(pa.pocketsinfo.GetMSENPR(i,80)):>10.4f}\n")

class MGHelper:
    @staticmethod
    def WritePocketPair(mg:ModelGroup,file_path:str,score_cutoff:float=0.0,dist_cutoff:float=2.0,**kw)->None:
        fdir=os.path.dirname(file_path)
        if not os.path.exists(fdir):
            print('Directory does not exist, create directory...',end='')
            os.makedirs(fdir)
            print('done',end='')
        pnum,opid,pcluster=mg.PocketMatch(score_cutoff,dist_cutoff)
        pcumnum=np.cumsum(pnum)#每个Modle的口袋数量
        cluster_num=int(np.max(pcluster))
        pcluster_s=np.split(pcluster,pcumnum)#每个口袋对应的簇ID
        opid_s=np.split(opid,pcumnum)#每个口袋的原始ID
        model_num=len(mg.pa_list)
        need_write_cluster_id=np.argwhere(np.array([np.sum([c in pcluster_s[m] for m in range(model_num)]) for c in range(cluster_num)])>1).flatten()#提取出需要处理的簇ID

        with open(file_path,'w') as f:
            out='#'
            for i in range(model_num):
                out+=f'model{i};'
            f.writelines(out[:-1]+'\n')

            pockpair_id=opid_s[0][np.argwhere(pcluster_s[0]==3).flatten()]
            for i in need_write_cluster_id:
                out=''
                for m in range(model_num):
                    pockpair_id=np.array([mg.pa_list[m].pocketsinfo.GetRank(o)+1 for o in opid_s[m][np.argwhere(pcluster_s[m]==i).flatten()]])
                    out+=np.array2string(pockpair_id,separator=',')[1:-1]+';'
                f.writelines(out[:-1]+'\n')

def str2bool(s:str)->bool:
    t=s.lower()
    if t!='true' and t!='false':
        print(f'wrong input \'{t}\'')
    return t=='true'

def ParserGeneral(md)->dict:
    out={'dist_cutoff':3.0,
         'is_use_score':True,
         'score_cutoff':10.0,
         'percentile':80}
    keys=md.keys()
    if 'dist_cutoff' in keys:
        out['dist_cutoff']=float(md['dist_cutoff'])
        print(f"Set dist_cutoff={float(md['dist_cutoff'])} done.")
    
    if 'is_use_score' in keys:
        out['is_use_score']=bool(md['is_use_score'])
        print(f"Set is_use_score={bool(md['is_use_score'])} done.")

    if 'score_cutoff' in keys:
        out['score_cutoff']=float(md['score_cutoff'])
        print(f"Set score_cutoff={float(md['score_cutoff'])} done.")

    if 'percentile' in keys:
        out['percentile']=int(md['percentile'])
        print(f"Set percentile={int(md['percentile'])} done.")
    
    return out

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

def GetPA(md,gparam)->PocketsAnalysis:
    if 'unpickle' in md.keys() and len(md['unpickle'])!=0:
        print('Start deserializing data...',end='')
        ppath=md['unpickle']
        if not os.path.exists(ppath):
            print(f"file dose not exist! ({ppath})")
        with open(ppath,'rb') as f:
            pa=pickle.load(f)
        print(f'{"done":>8s}.')
        return pa

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
    lig_cutoff=4.0 #默认值修改位置
    rec_mask=[]
    if 'lig_mask' in md.keys() and len(md['lig_mask'])!=0:
        print('Parse  lig mask...',end='')
        lig_mask=ParserMask(protein,md['lig_mask'])
        if 'lig_cutoff' in md:
            lig_cutoff=float(md['lig_cutoff'])
            print(f'  ->set lig_cutoff({lig_cutoff})->  ',end='')
        print(f'{"done":>8s}.')
    if 'rec_mask' in md.keys() and len(md['rec_mask'])!=0:
        print('Parse  rec mask...',end='')
        rec_mask=ParserMask(protein,md['rec_mask'])
        print(f'{"done":>8s}.')
    pa=PocketsAnalysis(protein,rec_mask=rec_mask,lig_mask=lig_mask)
    pa.SetLigCutoff(lig_cutoff)

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

    pa.SetDistCutoff(gparam['dist_cutoff'])
    pa.isUseScore(gparam['is_use_score'])
    pa.SetScoreCutoff(gparam['score_cutoff'])
    pa.SetPercentile(gparam['percentile'])

    return pa

def WriteFiles(pa:PocketsAnalysis,md)->None:
    if 'out' not in md.keys():
        print("The 'out' parameter was not detected and the data file could not be exported")
        return
    out_dir=md['out']
    if str2bool(md.get('out_summary','True')):
        print('Write sub-pocket info...',end='')
        PAHelper.WriteSummaryInfo(pa,os.path.join(out_dir,'pock_info'))
        print(f'{"done":>8s}.')
    if str2bool(md.get('out_best','True')):
        print('Write best conformation...',end='')
        PAHelper.WriteBestFrames(pa,os.path.join(out_dir,'pock_info/best_conf'),life_time_coutoff=int(md.get('best_lt_cutoff','10')),frame_shift=int(md.get('frame_start','0')),reserve_num=int(md.get('best_num','10')),is_use_scoring=str2bool(md.get('best_use_score','True')))
        print(f'{"done":>8s}.')
    if str2bool(md.get('out_pock_info','False')):
        print('Write pocket information...',end='')
        PAHelper.WritePocketInfos(pa,os.path.join(out_dir,'pock_info/'),ftype=md.get('pock_info_ftype','npy'))
        print(f'{"done":>8s}.')
    if str2bool(md.get('out_main_best','False')):
        print('Write main group...',end='')
        PAHelper.WriteMainBestFrames(pa,os.path.join(out_dir,'pock_info/main_group'),life_time_coutoff=int(md.get('main_lt_cutoff','10')),frame_shift=int(md.get('frame_start','0')),reserve_num=int(md.get('main_num','1')),is_use_scoring=str2bool(md.get('main_use_score','True')))
        print(f'{"done":>8s}.')
    if str2bool(md.get('out_coex','False')):
        print('Write coexist & correlation matrix...',end='')
        PAHelper.WritePocketAccor(pa,os.path.join(out_dir,'pock_info/'),score_cutoff=float(md.get('score_cutoff','10')),coex_sort_by=md.get('coex_sort_by','id'),corr_sort_by=md.get('corr_sort_by','id'))
        print(f'{"done":>8s}.')
    if str2bool(md.get('out_group','False')):
        print('Write group info...',end='')
        PAHelper.WriteGroupVtime(pa,os.path.join(out_dir,'pock_info/'))
        print(f'{"done":>8s}.')
    if str2bool(md.get('out_mtp','False')):
        print('Write main group transformation probability matrix...',end='')
        PAHelper.WtriteMainTransProb(pa,os.path.join(out_dir,'pock_info/'),frame_cutoff=float(md.get('main_mtp_cutoff','10')))
        print(f'{"done":>8s}.\n')

def Serializing(pa,md)->None:
    if 'pickle' in md.keys() and len(md['pickle'])!=0:
        print('Start serializing data...',end='')
        ppath=md['pickle']
        if not os.path.exists(os.path.dirname(ppath)):
            print(f"dir dose not exist! ({os.path.dirname(ppath)})")
        with open(ppath,'wb') as f:
            pickle.dump(pa,f)
        print(f'{"done":>8s}.')

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

parser.add_argument('-v,--version', action='version', version='%(prog)s 1.2.3')
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
parser.add_argument('--lig_cutoff',type=float,default=4.0,help='Alpha spheres with a distance from the ligand exceeding this cutoff value will be deleted,defaulting to 4.0.')

parser.add_argument('--box',type=str,action='append',nargs='+',default=[],metavar='x',
                    help='''R|This parameter specifies the position boxes.
If you are only interested in a certain region of the protein, you can specify one or more location boxes that wrap around that region. 
The analysis will be focused within the location boxes, which allows the program output to filter out data you are not interested in.
At the same time, this can speed up the running of the program and reduce the memory requirements.
NOTE: To accommodate continuous changes in protein conformation, the box centers are set to the coordinates of a single atom or to the geometric centers of several atoms, 
which allows the box coordinates to be synchronized with the conformational changes of the protein.
Format: --box atom1_id,atom2_id length(X-axis direction),width(Y-axis direction),height(Z-axis direction)
Example: --box 372 18,10,22 458,963 14,12,20''')

parser.add_argument('--dist_cutoff',type=float,default=3.0,metavar='float, default=3.0',help='Distance cutoff value when clustering subpockets')


parser.add_argument('--frame_start',type=int,default=0,metavar='integer, default=0',help='This parameter specifies which frame of the trajectory to start processing from, and frames before this parameter will be ignored')
parser.add_argument('--frame_stop',type=int,default=-1,metavar='integer, default=-1',help='This parameter specifies the frame to which the program ends processing, and the frames after this parameter will be ignored. The default value of this parameter is the last frame')
parser.add_argument('--frame_offset',type=int,default=1,metavar='integer, default=1',help='This parameter specifies the frame interval step size for processing trajectories.')

parser.add_argument('--out',type=str,metavar='out_dir',default='./', help='The folder where the analysis results are saved')
parser.add_argument('--out_summary',type=bool,metavar='bool, default=True',default=True,help='Output the sub pocket info')
parser.add_argument('--out_best',type=bool,metavar='bool, default=True',default=True,help='Output the optimal conformation, which means the overall score of the pocket is the highest')
parser.add_argument('--out_pock_info',type=bool,metavar='bool, default=False',default=False,help='Save the original information of the pocket in. npy format. For further analysis')
parser.add_argument('--out_main_best',type=bool,metavar='bool, default=False',default=False,help='Output the optimal conformation of each main_group')
parser.add_argument('--out_coex',type=bool,metavar='bool, default=False',default=False,help='Control whether to output sub pocket coexistence and correlation matrix')
parser.add_argument('--out_group',type=bool,metavar='bool, default=False',default=False,help='Control whether to output sub pocket correlation matrix')
parser.add_argument('--out_mtp',type=bool,metavar='bool, default=False',default=False,help='Control whether to output main group transformation probability matrix')

parser.add_argument('--is_use_score',type=bool,default=True,metavar='integer, default=10, unit: frame',help='Whether to use pocket scoring function when calculating the optimal conformation. If it is not, the pocket volume will be used to rate the pocket.')
parser.add_argument('--percentile',type=int,help='Characterize the value of hyperparameter n in the sub pocket fluctuation amplitude function. An integer with values ranging from 0 to 100, defaulting to 80')
parser.add_argument('--score_cutoff',type=int,help='Score below score_cutoff sub pockets will not be included in the pocket scoring calculation. The default is 10.0.')
parser.add_argument('--best_num',type=int,default=10,metavar='integer, default=10',help='Maximum number of best conformation retention')
parser.add_argument('--main_num',type=int,default=1,help='Maximum number of best main group conformation retention.')
parser.add_argument('--pock_info_ftype',type=str,choices=['npy','txt'],default='npy',help='The file format for pocket information, default to npy.')
parser.add_argument('--coex_sort_by',type=str,default='id',choices=['id','rank'],help='In what order are the rows and columns of the coex matrix arranged')
parser.add_argument('--corr_sort_by',type=str,default='id',choices=['id','rank'],help='In what order are the rows and columns of the corr matrix arranged')

parser.add_argument('--config',type=str,default='',help='Specify the path to the control file. When this parameter is specified, all input information will be read from the control file. Other command line input parameters will be ignored')

parser.add_argument('--pickle',type=str,help='Serialize and save the analysis results')
parser.add_argument('--unpickle',type=str,help='Read saved data')

#%% run
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

    gparam=ParserGeneral(config['GENERAL'])

    if config['GENERAL']['mode']=='single':
        print('single mode')
        pa=GetPA(config['MODEL0'],gparam)
        if not 'unpickle' in config['MODEL0'].keys():
            ts=int(config[f'MODEL0'].get('frame_start','0'))
            te=int(config[f'MODEL0'].get('frame_stop','-1'))
            to=int(config[f'MODEL0'].get('frame_offset','1'))
            pa.Analysis(ts,te,to)
        WriteFiles(pa,config['MODEL0'])
        Serializing(pa,config['MODEL0'])
    elif config['GENERAL']['mode']=='multi':
        print('multi mode')
        if not 'pock_pairs' in config['GENERAL']:
                print('error! pock_pairs not set')
                exit()
        models=ModelGroup()
        align_id=[]
        step=[]

        for i in range(100):
            if f'MODEL{i}' in config:
                print(f'MODEL{i}')
                # for k,v in config[f'MODEL{i}'].items():
                #     print(f'{k} : {v}')
                models.AddPA(GetPA(config[f'MODEL{i}'],gparam))
                if not 'unpickle' in config[f'MODEL{i}'].keys():
                    models.need_analysis.append(i)
                if not 'align_mask' in config[f'MODEL{i}']:
                    print('error! align_mask not set')
                    exit()
                align_id.append(np.array(ParserMask(models[i].rec,config[f'MODEL{i}']['align_mask'].strip(),mode='bb'),dtype=int))
                ts=int(config[f'MODEL{i}'].get('frame_start','0'))
                te=int(config[f'MODEL{i}'].get('frame_stop','-1'))
                to=int(config[f'MODEL{i}'].get('frame_offset','1'))
                step.append([ts,te,to])
        print('Process model...')
        models.Analysis(step)
        print('process done\n')
        print('Align model...')
        # for k in models.pa_list:
        #     for m in range(k.size()):
        #         k.snap_shots[m]._pocket_xyz=np.stack(k.snap_shots[m]._pocket_xyz,axis=0)
        models.AlignModel(align_id)
        print(f'align done\nMatch pocket id...',end='')
        m_score=0.0
        m_dist=2.0
        if 'match_score_cutoff' in config['GENERAL']:
            m_score=float(config['GENERAL']['match_score_cutoff'])
            print(f'\nset "match_score_cutoff" to {m_score}')
        if 'match_dist_cutoff' in config['GENERAL']:
            m_dist=float(config['GENERAL']['match_dist_cutoff'])
            print(f'set "match_dist_cutoff" to {m_dist}')
        MGHelper.WritePocketPair(models,config['GENERAL']['pock_pairs'],m_score,m_dist)
        print(f'{"done":>8s}.\n')

        for i in range(len(models.pa_list)):
            print(f'Write data of MODEL{i}')
            WriteFiles(models.pa_list[i],config[f'MODEL{i}'])
            Serializing(models.pa_list[i],config[f'MODEL{i}'])
            print('')
