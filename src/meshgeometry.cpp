/*
 * meshgeometry.cpp
 *
 *  Created on: Dec 12, 2012
 *      Author: Marco Agnese, Francesco Boffa, Marco Chiaramello
 */

#include "meshgeometry.h"

void mesh_geometry::read_nodes_list(string file_name){

  fstream file;
  string my_str;
  istringstream my_str_converted;

  node temp_node;
  int cont(1);
  double x_temp,y_temp;

  file.open(file_name.c_str());

  if(file.is_open()){
    while(getline(file,my_str)){
      my_str_converted.str(my_str);
      my_str_converted>>x_temp>>y_temp;
      temp_node.set_i(cont);
      temp_node.set_x(x_temp);
      temp_node.set_y(y_temp);
      my_nodes.add_node(temp_node);
      my_str_converted.clear();
      ++cont;
    }
  }
}

void mesh_geometry::create_lines(void){
  my_lines.resize(my_nodes.num_nodes());
  for(int i=0;i!=my_lines.num_lines()-1;++i){
    my_lines.set_p(i,i+1,i);
    my_lines.set_i(i+1+my_nodes.num_nodes(),i);
    my_lines.set_lt(i+1,i);
  }
  my_lines.set_p(my_lines.num_lines()-1,0,my_lines.num_lines()-1);
  my_lines.set_i(my_lines.num_lines()+my_nodes.num_nodes(),my_lines.num_lines()-1);
  my_lines.set_lt(my_lines.num_lines(),my_lines.num_lines()-1);
}

void mesh_geometry::write_msh(string file_name){

  ofstream file;
  stringstream my_stringstream;

  file.open(file_name.c_str());

  if(file.is_open()) {
    file<<"$MeshFormat\n2 0 8 \n$EndMeshFormat\n$Nodes\n";
    file<<my_nodes.num_nodes()<<endl;
    for(int i=0;i!=my_nodes.num_nodes();++i)
      file<<my_nodes.i(i)<<" "<<my_nodes.x(i)<<" "<<my_nodes.y(i)<<" 0"<<endl;
    file<<"$EndNodes\n$Elements"<<endl;
    file<<my_nodes.num_nodes()+my_lines.num_lines()<<endl;
    for(int j=0;j!=my_nodes.num_nodes();++j) file<<my_nodes.i(j)<<" 15 3 0 "<<my_nodes.i(j)<<" 0 "<<my_nodes.i(j)<<endl;
    for(int j=0;j!=my_lines.num_lines();++j) file<<my_lines.i(j)<<" 1 3 0 "<<my_lines.lt(j)<<" 0 "<<my_lines.p(j)[0]+1<<" "<<my_lines.p(j)[1]+1<<endl;
    file<<"$EndElements"<<endl;
    file<<"$GhostElements"<<endl;
    file<<my_elements.num_ghosts()<<endl;
    file<<"$EndGhostElements"<<endl;
    file.close();
  }
}

void mesh_geometry::write_geo(string file_name){

  ofstream file;
  stringstream my_stringstream;

  file.open(file_name.c_str());

  if(file.is_open()) {
    for(int i=0;i!=my_nodes.num_nodes();++i) file<<"Point("<<my_nodes.i(i)<<") = {"<<my_nodes.x(i)<<", "<<my_nodes.y(i)<<", 0, 1e+22};"<<endl;
    for(int j=0;j!=my_lines.num_lines();++j) file<<"Line("<<my_lines.i(j)<<") = {"<<my_lines.p(j)[0]+1<<", "<<my_lines.p(j)[1]+1<<"};"<<endl;
  }
}

void mesh_geometry::read_mesh_struct(string file_name,const int& lines_of_param){

  fstream file;
  string my_str;
  istringstream my_str_converted;
  int cont(0),temp(0);
  double x_temp, y_temp;

  file.open(file_name.c_str());
  structure_param.resize(static_cast<unsigned>(2*lines_of_param));

  if(file.is_open()){
    while(getline(file,my_str)){
    if (cont<lines_of_param){
      my_str_converted.str(my_str);
      my_str_converted>>structure_param[2*cont]>>structure_param[2*cont+1];
      temp+=structure_param[2*cont]*structure_param[2*cont+1];
    }
    if (cont==lines_of_param) my_nodes.resize(temp);
    if (cont>(lines_of_param-1)){
      my_str_converted.str(my_str);
      my_str_converted>>x_temp>>y_temp;
      my_nodes.set_i(cont-lines_of_param+1,cont-lines_of_param);
      my_nodes.set_x(x_temp,cont-lines_of_param);
      my_nodes.set_y(y_temp,cont-lines_of_param);
    }
    my_str_converted.clear();
    ++cont;
    }
  }
}

void mesh_geometry::create_elements(const int& tag_value,const loop& flag){

  vector<int> n_elements(0);
  n_elements.resize(structure_param.size()/2+1,0);


  if(flag==CLOSED){
    for(vector<int>::iterator iter=structure_param.begin();iter!=structure_param.end();iter+=2)
        n_elements[(iter-structure_param.begin())/2+1]=n_elements[(iter-structure_param.begin())/2]+(*(iter)-1)*(*(iter+1));
    my_elements.resize(*(n_elements.end()-1));
    for(int i=0;i!=static_cast<int>(n_elements.size()-1);++i){
      for(int j=n_elements[i];j!=n_elements[i+1];++j){
        my_elements.set_i(j+1,j);
        my_elements.set_v(((j/(structure_param[2*i]-1)+1)*structure_param[2*i]+j%(structure_param[2*i]-1))%my_nodes.num_nodes(),0,j);
        my_elements.set_v(((j/(structure_param[2*i]-1)+1)*structure_param[2*i]+j%(structure_param[2*i]-1)+1)%my_nodes.num_nodes(),1,j);
        my_elements.set_v(((j/(structure_param[2*i]-1))*structure_param[2*i]+j%(structure_param[2*i]-1)+1)%my_nodes.num_nodes(),2,j);
        my_elements.set_v(((j/(structure_param[2*i]-1))*structure_param[2*i]+j%(structure_param[2*i]-1))%my_nodes.num_nodes(),3,j);
        my_elements.set_tag(tag_value,j);
      }
    }
  }
  else{
    for(vector<int>::iterator iter=structure_param.begin();iter!=(structure_param.end()-2);iter+=2)
      n_elements[(iter-structure_param.begin())/2+1]=n_elements[(iter-structure_param.begin())/2]+(*(iter)-1)*(*(iter+1));
    *(n_elements.end()-1)=n_elements[(structure_param.size()-2)/2]+(*(structure_param.end()-2)-1)*(*(structure_param.end()-1)-1);

    my_elements.resize(*(n_elements.end()-1));
    for(int i=0;i!=static_cast<int>(n_elements.size()-1);++i){
      for(int j=n_elements[i];j!=n_elements[i+1];++j){
        my_elements.set_i(j+1,j);
        my_elements.set_v((j/(structure_param[2*i]-1)+1)*structure_param[2*i]+j%(structure_param[2*i]-1),0,j);
        my_elements.set_v((j/(structure_param[2*i]-1)+1)*structure_param[2*i]+j%(structure_param[2*i]-1)+1,1,j);
        my_elements.set_v((j/(structure_param[2*i]-1))*structure_param[2*i]+j%(structure_param[2*i]-1)+1,2,j);
        my_elements.set_v((j/(structure_param[2*i]-1))*structure_param[2*i]+j%(structure_param[2*i]-1),3,j);
        my_elements.set_tag(tag_value,j);
      }
    }
  }
}

void mesh_geometry::split_mesh(string file_name){

  vector<mesh> splitted_mesh(0);
  vector<int> tags;
  vector<int>::iterator pos;
  mesh empty_mesh;

  // splits elements according the tags
  for(int j=0;j!=my_elements.num_elements();++j){
    pos=find(tags.begin(),tags.end(),my_elements.tag(j));
    if(pos==tags.end()){
      tags.push_back(my_elements.tag(j));
      splitted_mesh.push_back(empty_mesh);
      pos=tags.end()-1;
    }
    (splitted_mesh[pos-tags.begin()].get_elements()).add_element(my_elements.ele(j));
  }

  // copy all the nodes list on each splitted mesh
  for(int k=0;k!=static_cast<int>(splitted_mesh.size());++k){
    for(int j=0;j!=my_nodes.num_nodes();++j)(splitted_mesh[k].get_nodes()).add_node(my_nodes.nod(j));
  }

  // splits lines according the tags
  for(int k=0;k!=static_cast<int>(splitted_mesh.size());++k){
    for(int j=0;j!=my_lines.num_lines();++j){
      for(int l=0;l!=(splitted_mesh[k].get_elements()).num_elements();++l){
        for(int i=0;i!=4;++i){
          if((my_lines.p(j)[0]==(((splitted_mesh[k].get_elements()).v(l))[i]))&&(my_lines.p(j)[1]==(((splitted_mesh[k].get_elements()).v(l))[(i+1)%4])))
            (splitted_mesh[k].get_lines()).add_line(my_lines.lin(j));
        }
      }
    }
  }

  // write splitted mesh
  string new_file_name;
  stringstream my_stringstream;
  string split_n;
  for(int k=0;k!=static_cast<int>(splitted_mesh.size());++k){
    new_file_name=file_name;
    new_file_name.resize(file_name.size()-4);
    my_stringstream<<k;
    split_n=my_stringstream.str();
    new_file_name=new_file_name+"_"+split_n+"_split.msh";
    splitted_mesh[k].write_file(new_file_name);
  }
}

