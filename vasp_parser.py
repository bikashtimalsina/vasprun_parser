import xml.etree.ElementTree as ET
class vasprun_parser:
    def __init__(self,path):
        self.path=path
    def get_root(self):
        tree=ET.parse(self.path)
        root=tree.getroot()
        return root
    def generator(self):
        gen_node=self.get_root().getchildren()[0]
        return [(gnode_item.attrib['name'],gnode_item.text) for gnode_item in gen_node.iter('i')]
    def incar(self):
        inc=self.get_root().getchildren()[1]
        return [(inc_i.attrib['name'],inc_i.text) for inc_i in inc.iter('i')]
    def monkhorst_pack(self):
        kpoints=self.get_root().getchildren()[2]
        return [(kpattr.attrib['name'],kpattr.text.split()) for kpattr in kpoints.getchildren()[0]]
    def kpoints_list(self):
        kpoints=self.get_root().getchildren()[2]
        return [kpoints.getchildren()[1].findall('v')[i].text.split() for i in range(len(kpoints.getchildren()[1].\
                                                                                         findall('v')))]
    def kpoints_weight(self):
        kpoints=self.get_root().getchildren()[2]
        return [kpoints.getchildren()[2].findall('v')[i].text.strip() for i in range(len(kpoints.getchildren()[2].\
                                                                                         findall('v')))]
    def parameters(self):
        mylist=[]
        param=self.get_root().getchildren()[3]
        iall=param.iter('i')
        vall=param.iter('v')
        for iinfo in iall:
            mylist.append((iinfo.attrib['name'],iinfo.text.strip()))
        for vinfo in vall:
            mylist.append((vinfo.attrib['name'],vinfo.text.split()))
        return mylist
    def atoms_info(self):
        atomic_information=[]
        atinfo=self.get_root().getchildren()[4]
        atomic_information.append(('No. of atoms',atinfo.getchildren()[0].text.strip()))
        atomic_information.append(('atom types', atinfo.getchildren()[1].text.strip()))
        ctags=atinfo.iter('c')
        for ctag in ctags:
            atomic_information.append(ctag.text.strip())
        return atomic_information
    def structure(self):
        cstruct=[];
        atom_struct=self.get_root().getchildren()[5]
        crystal=atom_struct.findall('crystal')
        elem_varrays=crystal[0].findall('varray')
        elemi=crystal[0].findall('i')
        for i in range(len(elem_varrays)):
            cstruct.append((elem_varrays[i].attrib['name'],\
            [elem_varrays[0].findall('v')[j].text.split() for j in range(len(elem_varrays[0].findall('v')))]))
        cstruct.append((elemi[0].attrib['name'],elemi[0].text.strip()))
        for data in atom_struct.findall('varray'):
            cstruct.append((data.attrib['name'],[data.findall('v')[j].text.split() \
            for j in range(len(data.findall('v')))]))
        return cstruct
    def calculation(self):
        calc_res=[]
        calc_forces=[]
        calc_positions=[]
        calc=self.get_root().getchildren()
        for i in range(6,len(calc)):
            tag_list=[]
            itags=calc[i].iter('i')
            timetags=calc[i].iter('time')
            cpos=[]
            cforces=[]
            for itag in itags:
                tag_list.append((itag.attrib['name'],itag.text.strip()))
            for timetag in timetags:
                tag_list.append((timetag.attrib['name'],timetag.text.split()))
            cal_struct=calc[i].findall('structure')
            cal_crystal=cal_struct[0].findall('crystal')
            struct_varray=cal_struct[0].findall('varray')
            force_varray=calc[i].findall('varray')
            for vs in struct_varray[0].iter('v'):
                cpos.append(vs.text.split())
            for vs in force_varray[0].iter('v'):
                cforces.append(vs.text.split())
            calc_positions.append(cpos)
            calc_forces.append(cforces)
            calc_res.append(tag_list)
        return calc_res,calc_positions,calc_forces