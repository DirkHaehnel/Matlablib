%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  calcdipolepdb2r.m %%
%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculates the dipole moment of a protein from its 3D structure as contained in a 
%Protein Data Bank .pdb file (http://www.pdb.org). In this sense,
%the dipole moment is the first moment of the charge distribution around the rotational centre of 
%drag, which is estimated as the geometric centre of all atoms weighted by their van der Waals radius.
%These radii and the pH dependent partial charges of all polarised atoms and ionised groups are taken 
%from the PARSE parameter set [1]. Missing hydrogen atoms are constructed according to the equivalent 
%bond lengths in small organic molecules [2]. The structure must be complete for this script to work: 
%i.e., all non-hydrogen atoms, should be present.
%
%[1] Sitkoff, Sharp and Honig: J. Phys. Chem. 98, 1978-1988 (1994).
%[2] Sutton: Tables of interatomic distances and configurations in molecules and ions. 
%            Chemical Society Special Publications (1958).
%
%Usage:   [dipolemoment, dipolemagnitude, netcharge]=calcdipolepdb2r(filename, includechain, pH)
%
%where:
%dipolemoment is the dipole moment vector in the coordinate frame of the PDB data, and 
%dipolemagnitude is the magnitude of the dipole moment, both in S.I. units: Cm. 
%                (To convert to Debyes, divide by 3.335664e-30.)
%netcharge is the total ionic charge.
%filename is the .pdb file containing the atomic coordinates.
%includechain is a string of characters, identifying the names of chains to include. 
%               e.g. 'A', 'B' or 'AB' according to the pdb file.
%pH is pH.
function [dipolemoment,dipolemagnitude,netcharge]=calcdipolepdb2r(filename,includechain,pH)

%function [dipolemoment,dipolemagnitude,netcharge]=calcdipolepdb2r(filename,includechain,pH,chain,hets,info)
% old command: [moment,coordinates,charge]=calcdipolepdb2(chain,hets,info,includechain,includehets,pH)
%first read pdb file:
[chain,hets,info]=readpdb2r(filename);

if info.chaintally=='0'
    includechain='0';
    disp('Null chain only!')
end
if ~isempty(find(ismember(includechain,info.chaintally)==0))
    disp('Error: includechain input contains chain identifiers not in PDB file')
    return
end

includehets=0;

% PARSE pKa values
pKa_TYR=9.98;
pKa_CYS=10.30;
pKa_LYS=10.60;
pKa_ARG=13.65;
pKa_HIS=6.95;
pKa_ASP=4.80;
pKa_GLU=4.88;
pKa_Nterm=pKa_LYS;
pKa_Cterm=pKa_ASP;

standardAAs={'ALA' 'ARG' 'ASN' 'ASP' 'CYS' 'GLN' 'GLU' 'GLY' 'HIS' 'ILE' 'LEU' 'LYS' ...
            'MET' 'PRO' 'PHE' 'SER' 'THR' 'TRP' 'TYR' 'VAL'};

charge=[]; %row 1 ---net zero,
           %row 2 ---difference when ionised, weighted by ionisation fraction as func of pH
coordinates=[]; %coordinates of all (partially) charged atoms
allcoordinates=[]; %coordinates of all atoms

%loop though selected chains
for i=find(ismember(info.chaintally,includechain)) %indices of includechains
    %startindex=length(charge(1,:)); ['chain: ' info.chaintally(i) ', startindex= 'num2str(startindex)]
    
    %1. N terminal hydrogens
    if ismember(chain(i).residue(1).name, standardAAs)    
            NtermHs=find(ismember({chain(i).residue(1).atom.name},{'1H' '2H' '3H' '1D' '2D' '3D'}));
            if chain(i).residue(1).name=='PRO'
                NtermHs=[NtermHs find(ismember({chain(1).residue(1).atom.name},{'CB'}))];
            end
            if length(NtermHs) < 3
                chain(i).residue(1).atom(end+1).name='cH'; %label constructed H
                chain(i).residue(1).atom(end).element='H';
                cHoccupancy=3-length(NtermHs);
                chain(i).residue(1).atom(end).occupancy=cHoccupancy;
                %chain(i).residue(1).atom(end).occupancy=cHoccupancy;
                chain(i).residue(1).atom(find(ismember({chain(i).residue(1).atom.name},{'CA'})));
                ratoma=chain(i).residue(1).atom(find(ismember({chain(i).residue(1).atom.name},{'CA'}))).r;
                ratomb=chain(i).residue(1).atom(find(ismember({chain(i).residue(1).atom.name},{'N'}))).r;
                rba=ratoma-ratomb;
                uba=rba./(ones(3,1)*sqrt(dot(rba,rba)));
                chain(i).residue(1).atom(end).r=ratomb-uba*1.09/2;
            end
            NtermHs=[NtermHs find(ismember({chain(i).residue(1).atom.name},{'cH'}))];
            
            ionfracp=(1+10^(pH-pKa_Nterm))^-1;
            charge=[charge [+0.40/3*[chain(i).residue(1).atom(NtermHs).occupancy]; ... %+0.40 for hydrogen(s)
                    +1*ionfracp/3*[chain(i).residue(1).atom(NtermHs).occupancy]]]; %+1 amongst hydrogen(s)
            coordinates=[coordinates [chain(i).residue(1).atom(NtermHs).r]];
            charge(:,end+1)=[-0.40; 0]; %-0.40 for N
            coordinates(:,end+1)=chain(i).residue(1).atom(find(ismember({chain(i).residue(1).atom.name},{'N'}))).r;
            
            %plot these
		%     plot3(coordinates(1,:),coordinates(2,:),coordinates(3,:),'+')
		%     hold on
		%     plot3(ratoma(1),ratoma(2),ratoma(3),'s')
		%     plot3(ratomb(1),ratomb(2),ratomb(3),'^')
		%     plot3(ratomb(1)-uba(1)*1.09/2,ratomb(2)-uba(2)*1.09/2,ratomb(3)-uba(3)*1.09/2,'+r')
    end
    
    %2. C terminal HXT
    if ismember(chain(i).residue(end).name, standardAAs)    
            CtermC=find(ismember({chain(i).residue(end).atom.name},{'C'}));
            CtermN=find(ismember({chain(i).residue(end).atom.name},{'N'}));
            CtermCA=find(ismember({chain(i).residue(end).atom.name},{'CA'}));
            if isempty(CtermC)+isempty(CtermN)+isempty(CtermCA) %ignore terminus if any required atoms are absent
                disp(['Warning: atom(s) missing at C terminus ' chain(i).residue(end).name]);
            else
            CtermHXT=find(ismember({chain(i).residue(end).atom.name},{'HXT'}));
            CtermOXT=find(ismember({chain(i).residue(end).atom.name},{'OXT'}));
            CtermO=find(ismember({chain(i).residue(end).atom.name},{'O'}));
            if isempty(CtermO)
                chain(i).residue(end).atom(end+1).name='cO'; %label constructed O...in NCA direction
                chain(i).residue(end).atom(end).element='O';
                chain(i).residue(end).atom(end).occupancy=1;
                ratoma=chain(i).residue(end).atom(CtermN).r;
                ratomb=chain(i).residue(end).atom(CtermCA).r;
                ratomc=chain(i).residue(end).atom(CtermC).r;
                rba=ratoma-ratomb;
                uba=rba./(ones(3,1)*sqrt(dot(rba,rba)));
                chain(i).residue(end).atom(end).r=ratomc-uba*1.31;
            end
            CtermO=[CtermO find(ismember({chain(i).residue(end).atom.name},{'cO'}))];
            if isempty(CtermOXT)
                chain(i).residue(end).atom(end+1).name='cOXT'; %label constructed OXT
                chain(i).residue(end).atom(end).element='O';
                chain(i).residue(end).atom(end).occupancy=1;
                ratomA=chain(i).residue(end).atom(CtermO).r;
                ratomB=chain(i).residue(end).atom(CtermCA).r;
                ratomC=chain(i).residue(end).atom(CtermC).r;
		
                rCA=ratomA-ratomC;
                uCA=rCA./(ones(3,1)*sqrt(dot(rCA,rCA)));
                rCB=ratomB-ratomC;
                uCB=rCB./(ones(3,1)*sqrt(dot(rCB,rCB)));
                P=-1.25/sqrt(2*(1+dot(uCA,uCB)))*[1;1];
                
                chain(i).residue(end).atom(end).r=ratomC+P(1)*uCA+P(2)*uCB;
            end
            CtermOXT=[CtermOXT find(ismember({chain(i).residue(end).atom.name},{'cOXT'}))];
            if isempty(CtermHXT)
                chain(i).residue(end).atom(end+1).name='cHXT'; %label constructed HXT
                chain(i).residue(end).atom(end).element='H';
                chain(i).residue(end).atom(end).occupancy=1;
                ratomA=chain(i).residue(end).atom(CtermC).r;
                ratomB=chain(i).residue(end).atom(CtermO).r;
                ratomC=chain(i).residue(end).atom(CtermOXT).r;
		
                rCA=ratomA-ratomC;
                uCA=rCA./(ones(3,1)*sqrt(dot(rCA,rCA)));
                rCB=ratomA-ratomB;
                uCB=rCB./(ones(3,1)*sqrt(dot(rCB,rCB)));
                P=-0.95/sqrt(2*(1+dot(uCA,uCB)))*[1;1]; 
                
                chain(i).residue(end).atom(end).r=ratomC+P(1)*uCA+P(2)*uCB;
            end
            CtermHXT=[CtermHXT find(ismember({chain(i).residue(end).atom.name},{'cHXT'}))];
            
            ionfracn=(10^(pKa_Cterm-pH)+1)^-1;
            charge(:,end+1)=[+0.435; -0.435*ionfracn];% for HXT copying ASP
            coordinates(:,end+1)=chain(i).residue(end).atom(CtermHXT).r;
            charge(:,end+1)=[-0.49; -0.06*ionfracn]; % for OXT copying ASP
            coordinates(:,end+1)=chain(i).residue(end).atom(CtermOXT).r;
            charge(:,end+1)=[-0.495; -0.055*ionfracn]; % for O copying ASP
            coordinates(:,end+1)=chain(i).residue(end).atom(CtermO).r;
            charge(:,end+1)=[+0.55; -0.45*ionfracn]; % for C copying ASP
            coordinates(:,end+1)=chain(i).residue(end).atom(CtermC).r;
		
            %plot these
		%     plot3(coordinates(1,end-3:end),coordinates(2,end-3:end),coordinates(3,end-3:end),'+')
		%     hold on
		%     plot3(ratomA(1),ratomA(2),ratomA(3),'s')
		%     plot3(ratomB(1),ratomB(2),ratomB(3),'^')
		%     plot3(ratomC(1),ratomC(2),ratomC(3),'o')
    end
    end
    
    %3. bb hydrogens
    for j=[2:length(chain(i).resseqtally)]
    if ismember(chain(i).residue(j).name, standardAAs)            
        %find H (or CD for proline)
        bbC=find(ismember({chain(i).residue(j).atom.name},{'C'}));
        bbCA=find(ismember({chain(i).residue(j).atom.name},{'CA'}));
        bbN=find(ismember({chain(i).residue(j).atom.name},{'N'}));
        if isempty(bbC)+isempty(bbCA)+isempty(bbN) %ignore bb repeat if any required atoms are absent
            disp(['Warning: backbone atom(s) missing at ' chain(i).residue(j).name ...
                     ' ' num2str(chain(i).residue(j).seq)]);
        else
        bbH=find(ismember({chain(i).residue(j).atom.name},{'H' 'D'}));
        bbO=find(ismember({chain(i).residue(j).atom.name},{'O' 'cO'}));
        if strcmp(chain(i).residue(j).name,'PRO')
            bbH=find(ismember({chain(i).residue(j).atom.name},{'CD'}));
        end
        if isempty(bbH)
            chain(i).residue(j).atom(end+1).name='cH'; %label constructed H
            chain(i).residue(j).atom(end).element='H';
            chain(i).residue(j).atom(end).occupancy=1;
            ratomA=chain(i).residue(j).atom(bbC(1)).r;
            ratomB=chain(i).residue(j).atom(bbCA(1)).r;
            ratomC=chain(i).residue(j).atom(bbN(1)).r;
            rCA=ratomA-ratomC;
            uCA=rCA./(ones(3,1)*sqrt(dot(rCA,rCA)));
            rCB=ratomB-ratomC;
            uCB=rCB./(ones(3,1)*sqrt(dot(rCB,rCB)));
            P=-1.09/sqrt(2*(1+dot(uCA,uCB)))*[1;1]; 
            chain(i).residue(j).atom(end).r=ratomC+P(1)*uCA+P(2)*uCB;
        end
        bbH=[bbH find(ismember({chain(i).residue(j).atom.name},{'cH'}))];
        if isempty(bbO)
            ['*** Warning: constructed backbone O at residue index ' num2str(j) ' ***']
            chain(i).residue(j).atom(end+1).name='cO'; %label constructed O
            chain(i).residue(j).atom(end).element='O';
            chain(i).residue(j).atom(end).occupancy=1;
            ratoma=chain(i).residue(j).atom(bbN(1)).r;
            ratomb=chain(i).residue(j).atom(bbCA(1)).r;
            ratomc=chain(i).residue(j).atom(bbC(1)).r;
            rba=ratoma-ratomb;
            uba=rba./(ones(3,1)*sqrt(dot(rba,rba)));
            chain(i).residue(j).atom(end).r=ratomc-uba*1.25;
        end
        bbO=find(ismember({chain(i).residue(j).atom.name},{'O' 'cO'}));
        charge(:,end+1)=[+0.40; 0]; %+0.40 for H
        coordinates(:,end+1)=chain(i).residue(j).atom(bbH(1)).r;
        charge(:,end+1)=[-0.40; 0]; %-0.40 for N
        coordinates(:,end+1)=chain(i).residue(j).atom(bbN(1)).r;
        charge(:,end+1)=[+0.55; 0]; %+0.55 for C
        coordinates(:,end+1)=chain(i).residue(j).atom(bbC(1)).r;
        charge(:,end+1)=[-0.55; 0]; %-0.55 for O
        coordinates(:,end+1)=chain(i).residue(j).atom(bbO(1)).r;
    end
    %plot these
%      hold on
%      length(coordinates(1,:))
%      %temp=length(chain(i).resseqtally);
%      temp=2
%      jj=[(length(coordinates(1,:))-4*temp+1)+4:4:length(coordinates(1,:))];
%      plot3(coordinates(1,jj),coordinates(2,jj),coordinates(3,jj),'.b') %H
%      plot3(coordinates(1,jj+1),coordinates(2,jj+1),coordinates(3,jj+1),'.r') %N
%      plot3(coordinates(1,jj+2),coordinates(2,jj+2),coordinates(3,jj+2),'.g') %C
%      plot3(coordinates(1,jj+3),coordinates(2,jj+3),coordinates(3,jj+3),'.m') %O
    end
    end
    
    %4. side chain hydrogens -- scan all residues, pick out MET TYR TRP SER THR CYS ASN GLN LYS ARG HIS ASP GLU
    
    for j=find(ismember({chain(i).residue.name},'MET'))
        CG=find(ismember({chain(i).residue(j).atom.name},{'CG'}));
        SD=find(ismember({chain(i).residue(j).atom.name},{'SD'}));
        CE=find(ismember({chain(i).residue(j).atom.name},{'CE'}));
        if isempty(CG)+isempty(SD)+isempty(CE) %ignore residue if any required atoms are absent
            disp(['Warning: ignoring incomplete residue sidechain ' chain(i).residue(j).name ' ' num2str(chain(i).residue(j).seq)]);
        else
        charge(:,end+1)=[+0.265;0]; %+0.265 for CG
        coordinates(:,end+1)=chain(i).residue(j).atom(CG).r;
        charge(:,end+1)=[-0.53;0]; %-0.53 for SD
        coordinates(:,end+1)=chain(i).residue(j).atom(SD).r;
        charge(:,end+1)=[+0.265;0]; %+0.265 for CE
        coordinates(:,end+1)=chain(i).residue(j).atom(CE).r;
        end
    end

    for j=find(ismember({chain(i).residue.name},'TYR'))
        CB=find(ismember({chain(i).residue(j).atom.name},{'CB'}));
        CG=find(ismember({chain(i).residue(j).atom.name},{'CG'}));
        CD1=find(ismember({chain(i).residue(j).atom.name},{'CD1'}));
        CD2=find(ismember({chain(i).residue(j).atom.name},{'CD2'}));
        CE1=find(ismember({chain(i).residue(j).atom.name},{'CE1'}));
        CE2=find(ismember({chain(i).residue(j).atom.name},{'CE2'}));
        CZ=find(ismember({chain(i).residue(j).atom.name},{'CZ'}));
        OH=find(ismember({chain(i).residue(j).atom.name},{'OH'}));
        if isempty(CB)+isempty(CG)+isempty(CD1)+isempty(CD2)+isempty(CE1)+isempty(CE2)+isempty(CZ)+isempty(OH) %ignore residue if any required atoms are absent
            disp(['Warning: ignoring incomplete residue sidechain ' chain(i).residue(j).name ' ' num2str(chain(i).residue(j).seq)]);
        else
        HH=find(ismember({chain(i).residue(j).atom.name},{'HH'}));
        if isempty(HH) %construct ABH type
            chain(i).residue(j).atom(end+1).name='cHH'; %label constructed HH
            chain(i).residue(j).atom(end).element='H';
            chain(i).residue(j).atom(end).occupancy=1;
            ratoma=chain(i).residue(j).atom(CZ).r;
            ratomb=chain(i).residue(j).atom(OH).r;
            rba=ratoma-ratomb;
            %uba=rba./(ones(3,1)*sqrt(dot(rba,rba)));
            uba=rba/sqrt(dot(rba,rba));
            chain(i).residue(j).atom(end).r=ratomb-uba*0.95/2;
        end
        HH=[HH find(ismember({chain(i).residue(j).atom.name},{'cHH'}))];
        
        ionfracn=(10^(pKa_TYR-pH)+1)^-1;
        charge(:,end+1)=[+0.125; 0]; %+0.125 for CB, ion 0
        coordinates(:,end+1)=chain(i).residue(j).atom(CB).r;
        charge(:,end+1)=[-0.125; -0.07*ionfracn]; %-0.125 for CG, ionisation -0.07
        coordinates(:,end+1)=chain(i).residue(j).atom(CG).r;
        charge(:,end+1)=[0; -0.07*ionfracn]; %0 for CD1, ionisation -0.07
        coordinates(:,end+1)=chain(i).residue(j).atom(CD1).r;
        charge(:,end+1)=[0; -0.07*ionfracn]; %0 for CD2, ionisation -0.07
        coordinates(:,end+1)=chain(i).residue(j).atom(CD2).r;
        charge(:,end+1)=[0; -0.07*ionfracn]; %0 for CE1, ionisation -0.07
        coordinates(:,end+1)=chain(i).residue(j).atom(CD1).r;
        charge(:,end+1)=[0; -0.07*ionfracn]; %0 for CE2, ionisation -0.07
        coordinates(:,end+1)=chain(i).residue(j).atom(CD2).r;
        charge(:,end+1)=[+0.055; -0.205*ionfracn]; %+0.055 for CZ, ion -0.205
        coordinates(:,end+1)=chain(i).residue(j).atom(CZ).r;
        charge(:,end+1)=[-0.49; -0.01*ionfracn]; %-0.49 for OH, ion -0.01
        coordinates(:,end+1)=chain(i).residue(j).atom(OH).r;
        charge(:,end+1)=[+0.435; -0.435*ionfracn]; %+0.435 for HH, ion -0.435
        coordinates(:,end+1)=chain(i).residue(j).atom(HH).r;
    %plot these
    % plot3(coordinates(1,end+[-4:-1]),coordinates(2,end+[-4:-1]),coordinates(3,end+[-4:-1]),'+')
    % hold on; plot3(coordinates(1,end),coordinates(2,end),coordinates(3,end),'.g')
    end
    end

    for j=find(ismember({chain(i).residue.name},'TRP'))
        CG=find(ismember({chain(i).residue(j).atom.name},{'CG'}));
        CD2=find(ismember({chain(i).residue(j).atom.name},{'CD2'}));
        NE1=find(ismember({chain(i).residue(j).atom.name},{'NE1'}));
        CE2=find(ismember({chain(i).residue(j).atom.name},{'CE2'}));
        CE3=find(ismember({chain(i).residue(j).atom.name},{'CE3'}));
        CZ2=find(ismember({chain(i).residue(j).atom.name},{'CZ2'}));
        CZ3=find(ismember({chain(i).residue(j).atom.name},{'CZ3'}));
        CH2=find(ismember({chain(i).residue(j).atom.name},{'CH2'}));
        if isempty(CG)+isempty(CD2)+isempty(NE1)+isempty(CE2)+isempty(CE3)+isempty(CZ2)+isempty(CZ3)+isempty(CH2) %ignore residue if any required atoms are absent
            disp(['Warning: ignoring incomplete residue sidechain ' chain(i).residue(j).name ' ' num2str(chain(i).residue(j).seq)]);
        else
        HE1=find(ismember({chain(i).residue(j).atom.name},{'HE1'}));
        HE3=find(ismember({chain(i).residue(j).atom.name},{'HE3'}));
        HZ2=find(ismember({chain(i).residue(j).atom.name},{'HZ2'}));
        HZ3=find(ismember({chain(i).residue(j).atom.name},{'HZ3'}));
        HH2=find(ismember({chain(i).residue(j).atom.name},{'HH2'}));
        if isempty(HE1) %construct ABCH type
            chain(i).residue(j).atom(end+1).name='cHE1'; %label constructed HE1
            chain(i).residue(j).atom(end).element='H';
            chain(i).residue(j).atom(end).occupancy=1;
            ratomA=chain(i).residue(j).atom(CD1).r;
            ratomB=chain(i).residue(j).atom(CE2).r;
            ratomC=chain(i).residue(j).atom(NE1).r;
            rCA=ratomA-ratomC;
            uCA=rCA./(ones(3,1)*sqrt(dot(rCA,rCA)));
            rCB=ratomB-ratomC;
            uCB=rCB./(ones(3,1)*sqrt(dot(rCB,rCB)));
            P=-1.09/sqrt(2*(1+dot(uCA,uCB)))*[1;1]; 
            chain(i).residue(j).atom(end).r=ratomC+P(1)*uCA+P(2)*uCB;
        end
        HE1=[HE1 find(ismember({chain(i).residue(j).atom.name},{'cHE1'}))];
        if isempty(HZ2) %construct ABCH type
            chain(i).residue(j).atom(end+1).name='cHZ2'; %label constructed HZ2
            chain(i).residue(j).atom(end).element='H';
            chain(i).residue(j).atom(end).occupancy=1;
            ratomA=chain(i).residue(j).atom(CE2).r;
            ratomB=chain(i).residue(j).atom(CH2).r;
            ratomC=chain(i).residue(j).atom(CZ2).r;
            rCA=ratomA-ratomC;
            uCA=rCA./(ones(3,1)*sqrt(dot(rCA,rCA)));
            rCB=ratomB-ratomC;
            uCB=rCB./(ones(3,1)*sqrt(dot(rCB,rCB)));
            P=-1.08/sqrt(2*(1+dot(uCA,uCB)))*[1;1]; 
            chain(i).residue(j).atom(end).r=ratomC+P(1)*uCA+P(2)*uCB;
        end
        HZ2=[HZ2 find(ismember({chain(i).residue(j).atom.name},{'cHZ2'}))];
        if isempty(HH2) %construct ABCH type
            chain(i).residue(j).atom(end+1).name='cHH2'; %label constructed HH2
            chain(i).residue(j).atom(end).element='H';
            chain(i).residue(j).atom(end).occupancy=1;
            ratomA=chain(i).residue(j).atom(CZ2).r;
            ratomB=chain(i).residue(j).atom(CZ3).r;
            ratomC=chain(i).residue(j).atom(CH2).r;
            rCA=ratomA-ratomC;
            uCA=rCA./(ones(3,1)*sqrt(dot(rCA,rCA)));
            rCB=ratomB-ratomC;
            uCB=rCB./(ones(3,1)*sqrt(dot(rCB,rCB)));
            P=-1.08/sqrt(2*(1+dot(uCA,uCB)))*[1;1]; 
            chain(i).residue(j).atom(end).r=ratomC+P(1)*uCA+P(2)*uCB;
        end
        HH2=[HH2 find(ismember({chain(i).residue(j).atom.name},{'cHH2'}))];
        if isempty(HZ3) %construct ABCH type
            chain(i).residue(j).atom(end+1).name='cHZ3'; %label constructed HZ3
            chain(i).residue(j).atom(end).element='H';
            chain(i).residue(j).atom(end).occupancy=1;
            ratomA=chain(i).residue(j).atom(CH2).r;
            ratomB=chain(i).residue(j).atom(CE3).r;
            ratomC=chain(i).residue(j).atom(CZ3).r;
            rCA=ratomA-ratomC;
            uCA=rCA./(ones(3,1)*sqrt(dot(rCA,rCA)));
            rCB=ratomB-ratomC;
            uCB=rCB./(ones(3,1)*sqrt(dot(rCB,rCB)));
            P=-1.08/sqrt(2*(1+dot(uCA,uCB)))*[1;1]; 
            chain(i).residue(j).atom(end).r=ratomC+P(1)*uCA+P(2)*uCB;
        end
        HZ3=[HZ3 find(ismember({chain(i).residue(j).atom.name},{'cHZ3'}))];
        if isempty(HE3) %construct ABCH type
            chain(i).residue(j).atom(end+1).name='cHE3'; %label constructed HE3
            chain(i).residue(j).atom(end).element='H';
            chain(i).residue(j).atom(end).occupancy=1;
            ratomA=chain(i).residue(j).atom(CZ3).r;
            ratomB=chain(i).residue(j).atom(CD2).r;
            ratomC=chain(i).residue(j).atom(CE3).r;
            rCA=ratomA-ratomC;
            uCA=rCA./(ones(3,1)*sqrt(dot(rCA,rCA)));
            rCB=ratomB-ratomC;
            uCB=rCB./(ones(3,1)*sqrt(dot(rCB,rCB)));
            P=-1.08/sqrt(2*(1+dot(uCA,uCB)))*[1;1]; 
            chain(i).residue(j).atom(end).r=ratomC+P(1)*uCA+P(2)*uCB;
        end
        HE3=[HE3 find(ismember({chain(i).residue(j).atom.name},{'cHE3'}))];

        charge(:,end+1)=[+0.125; 0]; %+0.125 for CG
        coordinates(:,end+1)=chain(i).residue(j).atom(CG).r;
        charge(:,end+1)=[-0.125; 0]; %-0.125 for CD2
        coordinates(:,end+1)=chain(i).residue(j).atom(CD2).r;
        charge(:,end+1)=[-0.40; 0]; %-0.40 for NE1
        coordinates(:,end+1)=chain(i).residue(j).atom(NE1).r;
        charge(:,end+1)=[+0.40; 0]; %-0.49 for HE1
        coordinates(:,end+1)=chain(i).residue(j).atom(HE1).r;
        charge(:,end+1)=[0; 0]; %0 for CE2
        coordinates(:,end+1)=chain(i).residue(j).atom(CE2).r;
        charge(:,end+1)=[+0.125; 0]; %+0.125 for HZ2
        coordinates(:,end+1)=chain(i).residue(j).atom(HZ2).r;
        charge(:,end+1)=[-0.125; 0]; %-0.125 for CZ2
        coordinates(:,end+1)=chain(i).residue(j).atom(CZ2).r;
        charge(:,end+1)=[+0.125; 0]; %+0.125 for HH2
        coordinates(:,end+1)=chain(i).residue(j).atom(HH2).r;
        charge(:,end+1)=[-0.125; 0]; %-0.125 for CH2
        coordinates(:,end+1)=chain(i).residue(j).atom(CH2).r;
        charge(:,end+1)=[+0.125; 0]; %+0.125 for HZ3
        coordinates(:,end+1)=chain(i).residue(j).atom(HZ3).r;
        charge(:,end+1)=[-0.125; 0]; %-0.125 for CZ3
        coordinates(:,end+1)=chain(i).residue(j).atom(CZ3).r;
        charge(:,end+1)=[+0.125; 0]; %+0.125 for HE3
        coordinates(:,end+1)=chain(i).residue(j).atom(HE3).r;
        charge(:,end+1)=[-0.125; 0]; %-0.125 for CE3
        coordinates(:,end+1)=chain(i).residue(j).atom(CE3).r;
    %plot these
%     plot3(coordinates(1,end+[-11:0]),coordinates(2,end+[-11:0]),coordinates(3,end+[-11:0]),'+'); hold on
%     plot3(coordinates(1,end-12),coordinates(2,end-12),coordinates(3,end-12),'.g')
    end
    end

    for j=find(ismember({chain(i).residue.name},'SER'))
        CB=find(ismember({chain(i).residue(j).atom.name},{'CB'}));
        OG=find(ismember({chain(i).residue(j).atom.name},{'OG'}));
        if isempty(CB)+isempty(OG) %ignore residue if any required atoms are absent
            disp(['Warning: ignoring incomplete residue sidechain ' chain(i).residue(j).name ' ' num2str(chain(i).residue(j).seq)]);
        else
        HG=find(ismember({chain(i).residue(j).atom.name},{'HG'}));
        if isempty(HG) %construct ABH type
            chain(i).residue(j).atom(end+1).name='cHG'; %label constructed HG
            chain(i).residue(j).atom(end).element='H';
            chain(i).residue(j).atom(end).occupancy=1;
            ratoma=chain(i).residue(j).atom(CB).r;
            ratomb=chain(i).residue(j).atom(OG).r;
            rba=ratoma-ratomb;
            %uba=rba./(ones(3,1)*sqrt(dot(rba,rba)));
            uba=rba/sqrt(dot(rba,rba));
            chain(i).residue(j).atom(end).r=ratomb-uba*0.95/2;
        end
        HG=[HG find(ismember({chain(i).residue(j).atom.name},{'cHG'}))];
        charge(:,end+1)=[-0.49; 0];%-0.49 for OG
        coordinates(:,end+1)=chain(i).residue(j).atom(OG).r;
        charge(:,end+1)=[+0.49; 0];%+0.49 for HG
        coordinates(:,end+1)=chain(i).residue(j).atom(HG).r;
    %plot these
    % plot3(coordinates(1,end-1),coordinates(2,end-1),coordinates(3,end-1),'+')
    % hold on; plot3(coordinates(1,end),coordinates(2,end),coordinates(3,end),'.g')        
    % plot3(chain(i).residue(j).atom(CB).r(1),chain(i).residue(j).atom(CB).r(2),chain(i).residue(j).atom(CB).r(3),'s')
    end
    end
    
    for j=find(ismember({chain(i).residue.name},'THR'))
        CB=find(ismember({chain(i).residue(j).atom.name},{'CB'}));
        OG1=find(ismember({chain(i).residue(j).atom.name},{'OG1'}));
        if isempty(CB)+isempty(OG1) %ignore residue if any required atoms are absent
            disp(['Warning: ignoring incomplete residue sidechain ' chain(i).residue(j).name ' ' num2str(chain(i).residue(j).seq)]);
        else
        HG1=find(ismember({chain(i).residue(j).atom.name},{'HG1'}));
        if isempty(HG1) %construct ABH type
            chain(i).residue(j).atom(end+1).name='cHG1'; %label constructed HG1
            chain(i).residue(j).atom(end).element='H';
            chain(i).residue(j).atom(end).occupancy=1;
            ratoma=chain(i).residue(j).atom(CB).r;
            ratomb=chain(i).residue(j).atom(OG1).r;
            rba=ratoma-ratomb;
            %uba=rba./(ones(3,1)*sqrt(dot(rba,rba)));
            uba=rba/sqrt(dot(rba,rba));
            chain(i).residue(j).atom(end).r=ratomb-uba*0.95/2;
        end
        HG1=[HG1 find(ismember({chain(i).residue(j).atom.name},{'cHG1'}))];
        charge(:,end+1)=[-0.49; 0]; %-0.49 for OG1
        coordinates(:,end+1)=chain(i).residue(j).atom(OG1).r;
        charge(:,end+1)=[+0.49; 0]; %+0.49 for HG1
        coordinates(:,end+1)=chain(i).residue(j).atom(HG1).r;
    %plot these
    % plot3(coordinates(1,end-1),coordinates(2,end-1),coordinates(3,end-1),'+')
    % hold on; plot3(coordinates(1,end),coordinates(2,end),coordinates(3,end),'.g')        
    % plot3(chain(i).residue(j).atom(CB).r(1),chain(i).residue(j).atom(CB).r(2),chain(i).residue(j).atom(CB).r(3),'s')
    end
    end
    
    for j=find(ismember({chain(i).residue.name},'CYS'))
        CB=find(ismember({chain(i).residue(j).atom.name},{'CB'}));
        SG=find(ismember({chain(i).residue(j).atom.name},{'SG'}));
        if isempty(CB)+isempty(SG) %ignore residue if any required atoms are absent
            disp(['Warning: ignoring incomplete residue sidechain ' chain(i).residue(j).name ' ' num2str(chain(i).residue(j).seq)]);
        else
        HG=find(ismember({chain(i).residue(j).atom.name},{'HG'}));
        if isempty(HG) %construct ABH type
            chain(i).residue(j).atom(end+1).name='cHG'; %label constructed SG
            chain(i).residue(j).atom(end).element='H';
            chain(i).residue(j).atom(end).occupancy=1;
            ratoma=chain(i).residue(j).atom(CB).r;
            ratomb=chain(i).residue(j).atom(SG).r;
            rba=ratoma-ratomb;
            %uba=rba./(ones(3,1)*sqrt(dot(rba,rba)));
            uba=rba/sqrt(dot(rba,rba));
            chain(i).residue(j).atom(end).r=ratomb-uba*1.329/2;
        end
        HG=[HG find(ismember({chain(i).residue(j).atom.name},{'cHG'}))];

        ionfracn=(10^(pKa_CYS-pH)+1)^-1;
        charge(:,end+1)=[0; -0.08*ionfracn]; %0 for CB, ion -0.08
        coordinates(:,end+1)=chain(i).residue(j).atom(CB).r;
        charge(:,end+1)=[-0.29; -0.63*ionfracn]; %-0.29 for SG, ion -0.63
        coordinates(:,end+1)=chain(i).residue(j).atom(SG).r;
        charge(:,end+1)=[+0.29; -0.29*ionfracn]; %+0.29 for HG, ion -0.29
        coordinates(:,end+1)=chain(i).residue(j).atom(HG).r;
    %plot these
    % plot3(coordinates(1,end-1),coordinates(2,end-1),coordinates(3,end-1),'+')
    % hold on; plot3(coordinates(1,end),coordinates(2,end),coordinates(3,end),'.g')        
    % plot3(chain(i).residue(j).atom(CB).r(1),chain(i).residue(j).atom(CB).r(2),chain(i).residue(j).atom(CB).r(3),'s')
    end
    end
    
    for j=find(ismember({chain(i).residue.name},'ASN'))
        CG=find(ismember({chain(i).residue(j).atom.name},{'CG'}));
        OD1=find(ismember({chain(i).residue(j).atom.name},{'OD1'}));
        ND2=find(ismember({chain(i).residue(j).atom.name},{'ND2'}));
        onetwoHD2=find(ismember({chain(i).residue(j).atom.name},{'1HD2' '2HD2'}));
        if isempty(CG)+isempty(OD1)+isempty(ND2) %ignore residue if any required atoms are absent
            disp(['Warning: ignoring incomplete residue sidechain ' chain(i).residue(j).name ' ' num2str(chain(i).residue(j).seq)]);
        else
            if length(onetwoHD2) < 2 %construct ABH type
                chain(i).residue(j).atom(end+1).name='conetwo2HD2'; %label constructed onetwoHD2
                chain(i).residue(j).atom(end).element='H';
                chain(i).residue(j).atom(end).occupancy=2-length(onetwoHD2);
                ratoma=chain(i).residue(j).atom(CG).r;
                ratomb=chain(i).residue(j).atom(ND2).r;
                rba=ratoma-ratomb;
                %uba=rba./(ones(3,1)*sqrt(dot(rba,rba)));
                uba=rba/sqrt(dot(rba,rba));
                chain(i).residue(j).atom(end).r=ratomb-uba*1.09/2;
            end
            onetwoHD2=[onetwoHD2 find(ismember({chain(i).residue(j).atom.name},{'conetwo2HD2'}))];
            
            charge(:,end+1)=[-0.55; 0]; %-0.55 for OD1
            coordinates(:,end+1)=chain(i).residue(j).atom(OD1).r;
            charge(:,end+1)=[+0.55; 0]; %+0.55 for CG
            coordinates(:,end+1)=chain(i).residue(j).atom(CG).r;
            charge(:,end+1)=[-0.78; 0]; %-0.78 for ND2
            coordinates(:,end+1)=chain(i).residue(j).atom(ND2).r;
            charge(:,end+[1:length(onetwoHD2)])=[+0.78/2*[chain(i).residue(j).atom(onetwoHD2).occupancy]; 0*[chain(i).residue(j).atom(onetwoHD2).occupancy]]; %+0.78 for onetwoHD2(s)
            coordinates(:,end+[1:length(onetwoHD2)])=[chain(i).residue(j).atom(onetwoHD2).r];
          %plot these
	%       plot3(coordinates(1,end+[-2:0]-length(onetwoHD2)),coordinates(2,end+[-2:0]-length(onetwoHD2)),coordinates(3,end+[-2:0]-length(onetwoHD2)),'+')
	%       hold on; plot3(coordinates(1,end+[-length(onetwoHD2):-1]+1),coordinates(2,end+[-length(onetwoHD2):-1]+1),coordinates(3,end+[-length(onetwoHD2):-1]+1),'.g')        
        end
    end
    
    for j=find(ismember({chain(i).residue.name},'GLN'))
        CD=find(ismember({chain(i).residue(j).atom.name},{'CD'}));
        OE1=find(ismember({chain(i).residue(j).atom.name},{'OE1'}));
        NE2=find(ismember({chain(i).residue(j).atom.name},{'NE2'}));
        if isempty(CD)+isempty(OE1)+isempty(NE2) %ignore residue if any required atoms are absent
            disp(['Warning: ignoring incomplete residue sidechain ' chain(i).residue(j).name ' ' num2str(chain(i).residue(j).seq)]);
        else
        onetwoHE2=find(ismember({chain(i).residue(j).atom.name},{'1HE2' '2HE2'}));
        if length(onetwoHE2) < 2 %construct ABH type
            chain(i).residue(j).atom(end+1).name='conetwo2HE2'; %label constructed onetwoHE2
            chain(i).residue(j).atom(end).element='H';
            chain(i).residue(j).atom(end).occupancy=2-length(onetwoHE2);
            ratoma=chain(i).residue(j).atom(CD).r;
            ratomb=chain(i).residue(j).atom(NE2).r;
            rba=ratoma-ratomb;
            %uba=rba./(ones(3,1)*sqrt(dot(rba,rba)));
            uba=rba/sqrt(dot(rba,rba));
            chain(i).residue(j).atom(end).r=ratomb-uba*1.09/2;
        end
        onetwoHE2=[onetwoHE2 find(ismember({chain(i).residue(j).atom.name},{'conetwo2HE2'}))];
        
        charge(:,end+1)=[-0.55; 0]; %-0.55 for OE1
        coordinates(:,end+1)=chain(i).residue(j).atom(OE1).r;
        charge(:,end+1)=[+0.55; 0]; %+0.55 for CD
        coordinates(:,end+1)=chain(i).residue(j).atom(CD).r;
        charge(:,end+1)=[-0.78; 0]; %-0.78 for NE2
        coordinates(:,end+1)=chain(i).residue(j).atom(NE2).r;
        charge(:,end+[1:length(onetwoHE2)])=[+0.78/2*[chain(i).residue(j).atom(onetwoHE2).occupancy]; 0*[chain(i).residue(j).atom(onetwoHE2).occupancy]]; %+0.78 for onetwoHE2(s)
        coordinates(:,end+[1:length(onetwoHD2)])=[chain(i).residue(j).atom(onetwoHE2).r];
      %plot these
      % plot3(coordinates(1,end+[-2:0]-length(onetwoHE2)),coordinates(2,end+[-2:0]-length(onetwoHE2)),coordinates(3,end+[-2:0]-length(onetwoHE2)),'+')
      % hold on; plot3(coordinates(1,end+[-length(onetwoHE2):-1]+1),coordinates(2,end+[-length(onetwoHE2):-1]+1),coordinates(3,end+[-length(onetwoHE2):-1]+1),'.g')        
      end
    end
    
    for j=find(ismember({chain(i).residue.name},'LYS'))
        CE=find(ismember({chain(i).residue(j).atom.name},{'CE'}));
        NZ=find(ismember({chain(i).residue(j).atom.name},{'NZ'}));
        if isempty(CE)+isempty(NZ) %ignore residue if any required atoms are absent
            disp(['Warning: ignoring incomplete residue sidechain ' chain(i).residue(j).name ' ' num2str(chain(i).residue(j).seq)]);
        else
        onetwothreeHZ=find(ismember({chain(i).residue(j).atom.name},{'1HZ' '2HZ' '3HZ'}));
        if length(onetwothreeHZ) < 3 %construct ABH type
            chain(i).residue(j).atom(end+1).name='conetwothreeHZ'; %label constructed onetwothreeHZ
            chain(i).residue(j).atom(end).element='H';
            chain(i).residue(j).atom(end).occupancy=3-length(onetwothreeHZ);
            ratoma=chain(i).residue(j).atom(CE).r;
            ratomb=chain(i).residue(j).atom(NZ).r;
            rba=ratoma-ratomb;
            %uba=rba./(ones(3,1)*sqrt(dot(rba,rba)));
            uba=rba/sqrt(dot(rba,rba));
            chain(i).residue(j).atom(end).r=ratomb-uba*1.09/2;
        end
        onetwothreeHZ=[onetwothreeHZ find(ismember({chain(i).residue(j).atom.name},{'conetwothreeHZ'}))];

        ionfracp=(1+10^(pH-pKa_LYS))^-1;
        charge(:,end+1)=[0; +0.33*ionfracp]; %0 for CE, ion +0.33
        coordinates(:,end+1)=chain(i).residue(j).atom(CE).r;
        charge(:,end+1)=[-0.78; +0.46*ionfracp]; %-0.78 for NZ, ion +0.46
        coordinates(:,end+1)=chain(i).residue(j).atom(NZ).r;
        charge(:,end+[1:length(onetwothreeHZ)])=[+0.78/3*[chain(i).residue(j).atom(onetwothreeHZ).occupancy];+0.21*ionfracp/3*[chain(i).residue(j).atom(onetwothreeHZ).occupancy]]; %+0.78 for onetwothreeHZ(s), ion +0.21
        coordinates(:,end+[1:length(onetwothreeHZ)])=[chain(i).residue(j).atom(onetwothreeHZ).r];
      %plot these
      % plot3(coordinates(1,end+[-1:0]-length(onetwothreeHZ)),coordinates(2,end+[-1:0]-length(onetwothreeHZ)),coordinates(3,end+[-1:0]-length(onetwothreeHZ)),'+')
      % hold on; plot3(coordinates(1,end+[-length(onetwothreeHZ):-1]+1),coordinates(2,end+[-length(onetwothreeHZ):-1]+1),coordinates(3,end+[-length(onetwothreeHZ):-1]+1),'.g')        
      end
    end
    
    for j=find(ismember({chain(i).residue.name},'ARG'))
        CD=find(ismember({chain(i).residue(j).atom.name},{'CD'}));
        NE=find(ismember({chain(i).residue(j).atom.name},{'NE'}));
        HE=find(ismember({chain(i).residue(j).atom.name},{'HE'}));
        CZ=find(ismember({chain(i).residue(j).atom.name},{'CZ'}));
        NH1=find(ismember({chain(i).residue(j).atom.name},{'NH1'}));
        NH2=find(ismember({chain(i).residue(j).atom.name},{'NH2'}));
        if isempty(CD)+isempty(NE)+isempty(CZ)+isempty(NH1)+isempty(NH2) %ignore residue if any required atoms are absent
            disp(['Warning: ignoring incomplete residue sidechain ' chain(i).residue(j).name ' ' num2str(chain(i).residue(j).seq)]);
        else
        onetwoHH1=find(ismember({chain(i).residue(j).atom.name},{'1HH1' '2HH1'}));
        onetwoHH2=find(ismember({chain(i).residue(j).atom.name},{'1HH2' '2HH2'}));
        if isempty(HE) %construct ABCH type
            chain(i).residue(j).atom(end+1).name='cHE'; %label constructed HE
            chain(i).residue(j).atom(end).element='H';
            chain(i).residue(j).atom(end).occupancy=1;
            ratomA=chain(i).residue(j).atom(CZ).r;
            ratomB=chain(i).residue(j).atom(CD).r;
            ratomC=chain(i).residue(j).atom(NE).r;
            rCA=ratomA-ratomC;
            uCA=rCA./(ones(3,1)*sqrt(dot(rCA,rCA)));
            rCB=ratomB-ratomC;
            uCB=rCB./(ones(3,1)*sqrt(dot(rCB,rCB)));
            P=-1.09/2/sqrt(2*(1+dot(uCA,uCB)))*[1;1]; 
            chain(i).residue(j).atom(end).r=ratomC+P(1)*uCA+P(2)*uCB;
        end
        HE=[HE find(ismember({chain(i).residue(j).atom.name},{'cHE'}))];
        if length(onetwoHH1) < 2 %construct ABH type
            chain(i).residue(j).atom(end+1).name='conetwoHH1'; %label constructed onetwoHH1
            chain(i).residue(j).atom(end).element='H';
            chain(i).residue(j).atom(end).occupancy=2-length(onetwoHH1);
            ratoma=chain(i).residue(j).atom(CZ).r;
            ratomb=chain(i).residue(j).atom(NH1).r;
            rba=ratoma-ratomb;
            %uba=rba./(ones(3,1)*sqrt(dot(rba,rba)));
            uba=rba/sqrt(dot(rba,rba));
            chain(i).residue(j).atom(end).r=ratomb-uba*1.09/2;
        end
        onetwoHH1=[onetwoHH1 find(ismember({chain(i).residue(j).atom.name},{'conetwoHH1'}))];
        if length(onetwoHH2) < 2 %construct ABH type
            chain(i).residue(j).atom(end+1).name='conetwoHH2'; %label constructed onetwoHH2
            chain(i).residue(j).atom(end).element='H';
            chain(i).residue(j).atom(end).occupancy=2-length(onetwoHH2);
            ratoma=chain(i).residue(j).atom(CZ).r;
            ratomb=chain(i).residue(j).atom(NH2).r;
            rba=ratoma-ratomb;
            %uba=rba./(ones(3,1)*sqrt(dot(rba,rba)));
            uba=rba/sqrt(dot(rba,rba));
            chain(i).residue(j).atom(end).r=ratomb-uba*1.09/2;
        end
        onetwoHH2=[onetwoHH2 find(ismember({chain(i).residue(j).atom.name},{'conetwoHH2'}))];
        
        ionfracp=(1+10^(pH-pKa_ARG))^-1;
        charge(:,end+1)=[+0.28; +0.07*ionfracp]; %+0.28 for CD, ion +0.07
        coordinates(:,end+1)=chain(i).residue(j).atom(CD).r;
        charge(:,end+1)=[-0.56; +0.21*ionfracp]; %-0.56 for NE, ion +0.21
        coordinates(:,end+1)=chain(i).residue(j).atom(NE).r;
        charge(:,end+1)=[0; +0.45*ionfracp]; %0 for HE, ion +0.45
        coordinates(:,end+1)=chain(i).residue(j).atom(HE).r;
        charge(:,end+1)=[+0.28; +0.07*ionfracp]; %+0.28 for CZ, ion +0.07
        coordinates(:,end+1)=chain(i).residue(j).atom(CZ).r;
        charge(:,end+1)=[-0.75; +0.05*ionfracp]; %-0.75 for NH1, ion +0.05
        coordinates(:,end+1)=chain(i).residue(j).atom(NH1).r;
        charge(:,end+1)=[-0.75; +0.05*ionfracp]; %-0.75 for NH2, ion +0.05
        coordinates(:,end+1)=chain(i).residue(j).atom(NH2).r;
        charge(:,end+[1:length(onetwoHH1)])=[+0.75/2*[chain(i).residue(j).atom(onetwoHH1).occupancy];+0.05*ionfracp/2*[chain(i).residue(j).atom(onetwoHH1).occupancy]]; %+0.75 for onetwoHH1(s), ion +0.05
        coordinates(:,end+[1:length(onetwoHH1)])=[chain(i).residue(j).atom(onetwoHH1).r];
        charge(:,end+[1:length(onetwoHH2)])=[+0.75/2*[chain(i).residue(j).atom(onetwoHH2).occupancy];+0.05*ionfracp/2*[chain(i).residue(j).atom(onetwoHH1).occupancy]]; %+0.75 for onetwoHH2(s), ion +0.05
        coordinates(:,end+[1:length(onetwoHH2)])=[chain(i).residue(j).atom(onetwoHH2).r];
      %plot these
      % plot3(coordinates(1,end+[-5:0]-length(onetwoHH1)-length(onetwoHH2)),coordinates(2,end+[-5:0]-length(onetwoHH1)-length(onetwoHH2)),coordinates(3,end+[-5:0]-length(onetwoHH1)-length(onetwoHH2)),'+')
      % hold on; plot3(coordinates(1,end+[(-length(onetwoHH1)-length(onetwoHH2)):-1]+1),coordinates(2,end+[(-length(onetwoHH1)-length(onetwoHH2)):-1]+1),coordinates(3,end+[(-length(onetwoHH1)-length(onetwoHH2)):-1]+1),'.g')
      end
    end
    
    for j=find(ismember({chain(i).residue.name},'HIS'))
        CG=find(ismember({chain(i).residue(j).atom.name},{'CG'}));
        ND1=find(ismember({chain(i).residue(j).atom.name},{'ND1'}));
        HD1=find(ismember({chain(i).residue(j).atom.name},{'HD1'}));
        CE1=find(ismember({chain(i).residue(j).atom.name},{'CE1'}));
        HE1=find(ismember({chain(i).residue(j).atom.name},{'HE1'}));
        NE2=find(ismember({chain(i).residue(j).atom.name},{'NE2'}));
        HE2=find(ismember({chain(i).residue(j).atom.name},{'HE2'}));
        CD2=find(ismember({chain(i).residue(j).atom.name},{'CD2'}));
        HD2=find(ismember({chain(i).residue(j).atom.name},{'HD2'}));
        if isempty(CG)+isempty(ND1)+isempty(CE1)+isempty(NE2)+isempty(CD2) %ignore residue if any required atoms are absent
            disp(['Warning: ignoring incomplete residue sidechain ' chain(i).residue(j).name ' ' num2str(chain(i).residue(j).seq)]);
        else

        if isempty(HD1) %construct ABCH type
            chain(i).residue(j).atom(end+1).name='cHD1'; %label constructed HD1
            chain(i).residue(j).atom(end).element='H';
            chain(i).residue(j).atom(end).occupancy=1;
            ratomA=chain(i).residue(j).atom(CG).r;
            ratomB=chain(i).residue(j).atom(CE1).r;
            ratomC=chain(i).residue(j).atom(ND1).r;
            rCA=ratomA-ratomC;
            uCA=rCA./(ones(3,1)*sqrt(dot(rCA,rCA)));
            rCB=ratomB-ratomC;
            uCB=rCB./(ones(3,1)*sqrt(dot(rCB,rCB)));
            P=-1.09/sqrt(2*(1+dot(uCA,uCB)))*[1;1]; 
            chain(i).residue(j).atom(end).r=ratomC+P(1)*uCA+P(2)*uCB;
        end
        HD1=[HD1 find(ismember({chain(i).residue(j).atom.name},{'cHD1'}))];
        if isempty(HE1) %construct ABCH type
            chain(i).residue(j).atom(end+1).name='cHE1'; %label constructed HE1
            chain(i).residue(j).atom(end).element='H';
            chain(i).residue(j).atom(end).occupancy=1;
            ratomA=chain(i).residue(j).atom(ND1).r;
            ratomB=chain(i).residue(j).atom(NE2).r;
            ratomC=chain(i).residue(j).atom(CE1).r;
            rCA=ratomA-ratomC;
            uCA=rCA./(ones(3,1)*sqrt(dot(rCA,rCA)));
            rCB=ratomB-ratomC;
            uCB=rCB./(ones(3,1)*sqrt(dot(rCB,rCB)));
            P=-1.08/sqrt(2*(1+dot(uCA,uCB)))*[1;1]; 
            chain(i).residue(j).atom(end).r=ratomC+P(1)*uCA+P(2)*uCB;
        end
        HE1=[HE1 find(ismember({chain(i).residue(j).atom.name},{'cHE1'}))];
        if isempty(HE2) %construct ABCH type
            chain(i).residue(j).atom(end+1).name='cHE2'; %label constructed HE2
            chain(i).residue(j).atom(end).element='H';
            chain(i).residue(j).atom(end).occupancy=1;
            ratomA=chain(i).residue(j).atom(CE1).r;
            ratomB=chain(i).residue(j).atom(CD2).r;
            ratomC=chain(i).residue(j).atom(NE2).r;
            rCA=ratomA-ratomC;
            uCA=rCA./(ones(3,1)*sqrt(dot(rCA,rCA)));
            rCB=ratomB-ratomC;
            uCB=rCB./(ones(3,1)*sqrt(dot(rCB,rCB)));
            P=-1.09/sqrt(2*(1+dot(uCA,uCB)))*[1;1]; 
            chain(i).residue(j).atom(end).r=ratomC+P(1)*uCA+P(2)*uCB;
        end
        HE2=[HE2 find(ismember({chain(i).residue(j).atom.name},{'cHE2'}))];
        if isempty(HD2) %construct ABCH type
            chain(i).residue(j).atom(end+1).name='cHD2'; %label constructed HD2
            chain(i).residue(j).atom(end).element='H';
            chain(i).residue(j).atom(end).occupancy=1;
            ratomA=chain(i).residue(j).atom(NE2).r;
            ratomB=chain(i).residue(j).atom(CG).r;
            ratomC=chain(i).residue(j).atom(CD2).r;
            rCA=ratomA-ratomC;
            uCA=rCA./(ones(3,1)*sqrt(dot(rCA,rCA)));
            rCB=ratomB-ratomC;
            uCB=rCB./(ones(3,1)*sqrt(dot(rCB,rCB)));
            P=-1.08/sqrt(2*(1+dot(uCA,uCB)))*[1;1]; 
            chain(i).residue(j).atom(end).r=ratomC+P(1)*uCA+P(2)*uCB;
        end
        HD2=[HD2 find(ismember({chain(i).residue(j).atom.name},{'cHD2'}))];
        
        ionfracp=(1+10^(pH-pKa_HIS))^-1;
        charge(:,end+1)=[-0.40; +0.05*ionfracp]; %-0.40 for ND1, ion +0.05
        coordinates(:,end+1)=chain(i).residue(j).atom(ND1).r;
        charge(:,end+1)=[+0.40; +0.05*ionfracp];%+0.40 for HD1, ion +0.05
        coordinates(:,end+1)=chain(i).residue(j).atom(HD1).r;
        charge(:,end+1)=[+0.155; +0.12*ionfracp];%+0.155 for CE1, ion +0.12
        coordinates(:,end+1)=chain(i).residue(j).atom(CE1).r;
        charge(:,end+1)=[+0.125; 0];%+0.125 for HE1, ion 0
        coordinates(:,end+1)=chain(i).residue(j).atom(HE1).r;
        charge(:,end+1)=[-0.56; +0.21*ionfracp];%-0.56 for NE2, ion +0.21
        coordinates(:,end+1)=chain(i).residue(j).atom(NE2).r;
        charge(:,end+1)=[0; +0.45*ionfracp];%0 for HE2, ion +0.45
        coordinates(:,end+1)=chain(i).residue(j).atom(HE2).r;
        charge(:,end+1)=[+0.155; +0.12*ionfracp];%+0.155 for CD2, ion +0.12
        coordinates(:,end+1)=chain(i).residue(j).atom(CD2).r;
        charge(:,end+1)=[+0.125; 0];%+0.125 for HD2, ion 0
        coordinates(:,end+1)=chain(i).residue(j).atom(HD2).r;
        end
    end
      %plot these
      % plot3(coordinates(1,end+[-7:0]),coordinates(2,end+[-7:0]),coordinates(3,end+[-7:0]),'+')
      % hold on; plot3(coordinates(1,end+[(-length(onetwoHH1)-length(onetwoHH2)):-1]+1),coordinates(2,end+[(-length(onetwoHH1)-length(onetwoHH2)):-1]+1),coordinates(3,end+[(-length(onetwoHH1)-length(onetwoHH2)):-1]+1),'.g')        
    
    for j=find(ismember({chain(i).residue.name},'ASP'))
        CG=find(ismember({chain(i).residue(j).atom.name},{'CG'}));
        OD1=find(ismember({chain(i).residue(j).atom.name},{'OD1'}));
        OD2=find(ismember({chain(i).residue(j).atom.name},{'OD2'}));
        if isempty(CG)+isempty(OD1)+isempty(OD2) %ignore residue if any required atoms are absent
            disp(['Warning: ignoring incomplete residue sidechain ' chain(i).residue(j).name ' ' num2str(chain(i).residue(j).seq)]);
        else

        HD2=find(ismember({chain(i).residue(j).atom.name},{'HD2'}));
        if isempty(HD2)
            chain(i).residue(j).atom(end+1).name='cHD2'; %label constructed HD2
            chain(i).residue(j).atom(end).element='H';
            chain(i).residue(j).atom(end).occupancy=1;
            ratomA=chain(i).residue(j).atom(OD1).r;
            ratomB=chain(i).residue(j).atom(CG).r;
            ratomC=chain(i).residue(j).atom(OD2).r;
	
            rCA=ratomA-ratomC;
            uCA=rCA./(ones(3,1)*sqrt(dot(rCA,rCA)));
            rCB=ratomB-ratomA;
            uCB=rCB./(ones(3,1)*sqrt(dot(rCB,rCB)));
            P=-0.95/sqrt(2*(1+dot(uCA,uCB)))*[1;1]; 
            
            chain(i).residue(j).atom(end).r=ratomC+P(1)*uCA+P(2)*uCB;
        end
        HD2=[HD2 find(ismember({chain(i).residue(j).atom.name},{'cHD2'}))];

        ionfracn=(10^(pKa_ASP-pH)+1)^-1;
        charge(:,end+1)=[+0.55; -0.45*ionfracn]; %+0.55 for CG, ion -0.45
        coordinates(:,end+1)=chain(i).residue(j).atom(CG).r;
        charge(:,end+1)=[-0.495; -0.055*ionfracn]; %-0.495 for OD1, ion -0.055
        coordinates(:,end+1)=chain(i).residue(j).atom(OD1).r;
        charge(:,end+1)=[-0.49; -0.06*ionfracn]; %-0.49 for OD2, ion -0.06
        coordinates(:,end+1)=chain(i).residue(j).atom(OD2).r;
        charge(:,end+1)=[+0.435; -0.435*ionfracn]; %+0.435 for HD2, ion -0.435
        coordinates(:,end+1)=chain(i).residue(j).atom(HD2).r;
        end
    end
    %plot these
    % plot3(coordinates(1,end+[-3:-1]),coordinates(2,end+[-3:-1]),coordinates(3,end+[-3:-1]),'+')
    % hold on; plot3(coordinates(1,end),coordinates(2,end),coordinates(3,end),'.g')
    
    for j=find(ismember({chain(i).residue.name},'GLU'))
        CD=find(ismember({chain(i).residue(j).atom.name},{'CD'}));
        OE1=find(ismember({chain(i).residue(j).atom.name},{'OE1'}));
        OE2=find(ismember({chain(i).residue(j).atom.name},{'OE2'}));
        if isempty(CD)+isempty(OE1)+isempty(OE2) %ignore residue if any required atoms are absent
            disp(['Warning: ignoring incomplete residue sidechain ' chain(i).residue(j).name ' ' num2str(chain(i).residue(j).seq)]);
        else
        HE2=find(ismember({chain(i).residue(j).atom.name},{'HE2'}));
        if isempty(HE2)
            chain(i).residue(j).atom(end+1).name='cHE2'; %label constructed HE2
            chain(i).residue(j).atom(end).element='H';
            chain(i).residue(j).atom(end).occupancy=1;
            ratomA=chain(i).residue(j).atom(OE1).r;
            ratomB=chain(i).residue(j).atom(CD).r;
            ratomC=chain(i).residue(j).atom(OE2).r;
	
            rCA=ratomA-ratomC;
            uCA=rCA./(ones(3,1)*sqrt(dot(rCA,rCA)));
            rCB=ratomB-ratomA;
            uCB=rCB./(ones(3,1)*sqrt(dot(rCB,rCB)));
            P=-0.95/sqrt(2*(1+dot(uCA,uCB)))*[1;1]; 
            
            chain(i).residue(j).atom(end).r=ratomC+P(1)*uCA+P(2)*uCB;
        end
        HE2=[HE2 find(ismember({chain(i).residue(j).atom.name},{'cHE2'}))];

        ionfracn=(10^(pKa_GLU-pH)+1)^-1;
        charge(:,end+1)=[+0.55; -0.45*ionfracn]; %+0.55 for CD, ion -0.45
        coordinates(:,end+1)=chain(i).residue(j).atom(CD).r;
        charge(:,end+1)=[-0.495; -0.055*ionfracn]; %-0.495 for OE1, ion -0.055
        coordinates(:,end+1)=chain(i).residue(j).atom(OE1).r;
        charge(:,end+1)=[-0.49; -0.06*ionfracn]; %-0.49 for OE2, ion -0.06
        coordinates(:,end+1)=chain(i).residue(j).atom(OE2).r;
        charge(:,end+1)=[+0.435; -0.435*ionfracn]; %+0.435 for HE2, ion -0.435
        coordinates(:,end+1)=chain(i).residue(j).atom(HE2).r;
        
        %plot these
        % plot3(coordinates(1,end+[-3:-1]),coordinates(2,end+[-3:-1]),coordinates(3,end+[-3:-1]),'+')
        % hold on; plot3(coordinates(1,end),coordinates(2,end),coordinates(3,end),'.g')
        end
        end
        
        %generate coordinate array of all atoms: vdW radius, occupancy, and coordinates
        chainatoms=[chain(i).residue.atom];
%        chainatomselement=[chainatoms.element]; % length(chainatomselement)
        chainatomsnames=char({chainatoms.name});
        chainatomselement=chainatomsnames(:,1);
        for istepalong=[2:size(chainatomsnames,2)]
            stepalong=~ismember(chainatomselement,'DHCNOS');
            chainatomselement(stepalong)=chainatomsnames(stepalong,istepalong);
        end
        if ~isempty(find(stepalong))
            disp('Error: atoms of indices below have unidentified element, how can that be?')
            disp(find(stepalong))
        end
        
        clear vdWradii
        vdWradii(find(chainatomselement=='N'))=1.5; % all according to PARSE
        vdWradii(find(chainatomselement=='O'))=1.4;
        vdWradii(find(chainatomselement=='C'))=1.7;
        vdWradii(find(chainatomselement=='S'))=1.85;
        vdWradii(find(chainatomselement=='H'))=1.0;
        vdWradii(find(chainatomselement=='D'))=1.0;
%         disp(i)
%         disp(size(vdWradii))
%         disp(size([chainatoms.occupancy;chainatoms.r]))
%         disp(size(allcoordinates))
%         disp(size([vdWradii; [chainatoms.occupancy;chainatoms.r]]))
        allcoordinates=[vdWradii; [chainatoms.occupancy;chainatoms.r]];
%        allcoordinates=[allcoordinates [vdWradii; [chainatoms.occupancy;chainatoms.r]]];
%         disp(size(allcoordinates))
end

% dipole moment of neutral component
p=sum(coordinates.*(ones(3,1)*charge(1,:)),2) ;
mag_p=sqrt(dot(p,p));
map_p_D=mag_p*1.6e-19*1e-10/3.335664e-30 ;

% centre of charge
centre_q=sum(coordinates.*(ones(3,1)*charge(2,:)),2)/sum(charge(2,:))*ones(1,length(coordinates(1,:)));
% for geometric centre, use allcoordinates coordinate array of all atoms and their vdW radii
allnatoms=length(allcoordinates(1,:));
% 1) geometric centre of all atoms, weighted by occupancy
centre_g1=sum(allcoordinates(3:5,:).*(ones(3,1)*allcoordinates(2,:)),2)/sum(allcoordinates(2,:))*ones(1,length(coordinates(1,:)));
% 2) vdW radius weighted geometric centre of all atoms, weighted by occupancy
centre_g2=sum(allcoordinates(3:5,:).*(ones(3,1)*(allcoordinates(1,:).*allcoordinates(2,:))),2)/sum(allcoordinates(1,:).*allcoordinates(2,:))*ones(1,length(coordinates(1,:)));

% 1) first moment (dipole moment) of charged component around simple geometric centre.
p_centre_g1 = (p + sum((coordinates-centre_g1).*(ones(3,1)*charge(2,:)), 2) )*1.6e-19*1e-10;
mag_p_centre_g1=sqrt(dot(p_centre_g1,p_centre_g1));
mag_p_centre_g1_D=mag_p_centre_g1/3.335664e-30 ;
% 2) first moment (dipole moment) of charged component around vdW radius weighted geometric centre.
p_centre_g2 = (p + sum((coordinates-centre_g2).*(ones(3,1)*charge(2,:)), 2) )*1.6e-19*1e-10;
mag_p_centre_g2=sqrt(dot(p_centre_g2,p_centre_g2));
mag_p_centre_g2_D=mag_p_centre_g2/3.335664e-30 ;

dipolemoment=p_centre_g2;
dipolemagnitude=mag_p_centre_g2 ;
netcharge=sum(charge(2,:));






function [chain,hets,info]=readpdb2r(filename)
%read entries in a pdb file:
%------------------------------
%read line by line, and take note of the following entries:
%
% HEADER    LIPOCALIN                               20-DEC-96   1BEB              
% TITLE     BOVINE BETA-LACTOGLOBULIN, LATTICE X                                  
% COMPND    MOL_ID: 1;                                                            
% COMPND   2 MOLECULE: BETA-LACTOGLOBULIN;                                        
% COMPND   3 CHAIN: A, B;                                                         
% COMPND   4 BIOLOGICAL_UNIT: PREDOMINANTLY DIMERIC                               
% ...
% ATOM      1  N   GLN A   5     -25.924  14.212  -4.336  1.00 27.22           N  
% ATOM      2  CA  GLN A   5     -26.461  15.615  -4.427  1.00 29.51           C  
% ...
% TER    1495      CYS A 160                                                      
% ...
% HETATM 2991  S   SO4 B   1      10.085  -3.150  17.197  0.64 63.26           S  
% HETATM 3000 1H   HOH   402      -1.277  26.680  -9.385  1.00 20.00           H  

fid=fopen(filename);
a=0;
info.compnd=[];
natoms=0; ii=0; ij=0; ik=0;
info.chaintally=[];
hets.hetseqtally=[];
j=0;
nhetatoms=0; ki=0; kj=0; kk=0;
info.header=[];
info.title=[];
info.compnd=[];
info.model=[];
while a~=1
    tline=fgetl(fid); 
    %disp(tline(1:5))
    findstr11toend=findstr(tline(11:end),'  ');
    if strncmp(tline,'HEADER ',7)
        info.header=[info.header tline(10+[1:findstr11toend(1)])];
    elseif strncmp(tline,'TITLE ',6)
        info.title=[info.title tline(10+[1:findstr11toend(1)])];
    elseif strncmp(tline,'COMPND ',7)
        info.compnd=[info.compnd tline(10+[1:findstr11toend(1)])];
    elseif strncmp(tline,'MODEL ',6)
        info.model=[info.model sscanf(tline(11:14),'%i')];
    elseif strncmp(tline,'ATOM ',5)*~ismember(tline(17),'BCDEFGHIJKLMNOPQRSTUVWXYZ')*~ismember(tline(27),'ABCDEFGHIJKLMNOPQRSTUVWXYZ')
        natoms=natoms+1;
        %
        chainidentifier=tline(22);
        if chainidentifier==' '
            chainidentifier='0'; %name a null chain
        end
        info.chaintally=unique([info.chaintally chainidentifier]); %tally of chains ABC...
        iiold=ii;
        ii=find(info.chaintally==chainidentifier); %index of chain
        chain(ii).name=chainidentifier;
        %
        resseq=sscanf(tline(23:26),'%i');
        resname=sscanf(tline(18:20),'%s');
        if(iiold<ii)
            chain(ii).resseqtally=[];
            chain(ii).residue(1).serialtally=[];
        end
        chain(ii).resseqtally=unique([chain(ii).resseqtally resseq]); %tally of res seqence numbers 5,6,7...
        ijold=ij;
        ij=find(chain(ii).resseqtally==resseq);   %index of residue in chain
        chain(ii).residue(ij).seq=resseq;
        chain(ii).residue(ij).name=resname;
        %
        serial=sscanf(tline(7:11),'%i');
        if(ijold<ij)
            chain(ii).residue(ij).serialtally=[];
        end
        chain(ii).residue(ij).serialtally=[chain(ii).residue(ij).serialtally serial];
        ik=find(chain(ii).residue(ij).serialtally==serial); %index of atom in residue in chain
        chain(ii).residue(ij).atom(ik).serial=serial;
        chain(ii).residue(ij).atom(ik).name=sscanf(tline(13:16),'%s');
        chain(ii).residue(ij).atom(ik).element=sscanf(tline(77:78),'%s');
        chain(ii).residue(ij).atom(ik).r=[sscanf(tline(31:38),'%e');sscanf(tline(39:46),'%e');sscanf(tline(47:54),'%e')];
        chain(ii).residue(ij).atom(ik).occupancy=sscanf(tline(55:60),'%e');
    elseif strncmp(tline,'TER ',4)
        j=j+1;
        ji=find(info.chaintally==chainidentifier);
        chain(ji).terserial=sscanf(tline(7:11),'%i');
        chain(ji).terseq=sscanf(tline(23:26),'%i');
    elseif strncmp(tline,'HETATM ',7)
        nhetatoms=nhetatoms+1;
        %
        ki=find(info.chaintally==chainidentifier); %index of chain if present
        %
        hetseq=sscanf(tline(23:26),'%i');
        hetname=sscanf(tline(18:20),'%s');
        hets.hetseqtally=unique([hets.hetseqtally hetseq]); %tally of het sequence numbers
        kjold=kj;
        kj=find(hets.hetseqtally==hetseq);   %index of het in sequence
        hets.het(kj).seq=hetseq;
        hets.het(kj).name=resname;
        hets.het(kj).chain=info.chaintally(ki);
        %
        serial=sscanf(tline(7:11),'%i');
        if(kjold<kj)
            hets.het(kj).serialtally=[];
        end
        hets.het(kj).serialtally=[hets.het(kj).serialtally serial];
        kk=find(hets.het(kj).serialtally==serial); %index of hetatom in het
        hets.het(kj).atom(kk).serial=serial;
        hets.het(kj).atom(kk).name=sscanf(tline(13:16),'%s');
        hets.het(kj).atom(kk).element=sscanf(tline(77:78),'%s');
        hets.het(kj).atom(kk).r=[sscanf(tline(31:38),'%e');sscanf(tline(39:46),'%e');sscanf(tline(47:54),'%e')];
        hets.het(kj).atom(kk).occupancy=sscanf(tline(55:60),'%e');
    elseif strncmp(tline,'END',3)
        a=1;
    elseif strncmp(tline,'ENDMDL',6)
        a=1;
    end
end
fclose(fid);

disp(['Title: ' info.title])
disp(['Header: ' info.header])
disp(['Compnd info: ' info.compnd])

