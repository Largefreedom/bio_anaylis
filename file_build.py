import os
import json

#打开json文件夹；
file = open('C:/Users/FREEDOM/Desktop/TCGA_data/mainfest/metadata.cart.2019-04-19.json', encoding='utf-8')
#读取json文件夹；
json_precess = json.loads(str(file.read()))
#创建空字典；
dict ={}
for i in json_precess:
    print(i['file_name'])
    print(i['associated_entities'][0]['entity_submitter_id'])
    dict[str(i['file_name']).strip('.gz')] =i['associated_entities'][0]['entity_submitter_id']
print(dict)


#mainfest文件路径；
path ='C:/Users/FREEDOM/Desktop/TCGA_data/mainfest'
filelist = os.listdir(path)#mainfest子路径下的所有文件列表；
for file_one in filelist:
    file = path + '/' +file_one
    print(file)
    list = os.listdir(file)[0]
    print(list)
    if '.gz' in list:
        olddir = file +'/' +list#原来文件名
        newdir = file + '/' + dict[list.split('.gz')[0]] +'.gz'#新的文件名
        os.rename(olddir,newdir)#重新命名
    # olddir = os.path.join(path,file_one)
    # print(olddir)
    # filename =os.path.splitext(file_one)[0]
    # newdir = os.path.join(path,dict[file_one])
    # os.rename(olddir,newdir)名

