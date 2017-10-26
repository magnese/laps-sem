/*
 * meshparallel.cpp
 *
 *  Created on: Jan 24, 2013
 *      Author: Marco Agnese, Francesco Boffa, Marco Chiaramello
 */

#include "meshparallel.h"

string mesh_parallel::file_name_pre(string file_name){

  stringstream my_stringstream;
  string rank_n;

  file_name.resize(file_name.size()-4);
  my_stringstream<<my_rank;
  rank_n=my_stringstream.str();
  file_name=file_name+"_"+rank_n+"_pre.msh";

  return file_name;
}

string mesh_parallel::file_name_par(string file_name){

  stringstream my_stringstream;
  string rank_n;

  file_name.resize(file_name.size()-4);
  my_stringstream<<my_rank;
  rank_n=my_stringstream.str();
  file_name=file_name+"_"+rank_n+"_par.msh";

  return file_name;
}

string mesh_parallel::file_name_ref(string file_name){

    stringstream my_stringstream;
    string rank_n;

    file_name.resize(file_name.size()-4);
    my_stringstream<<my_rank;
    rank_n=my_stringstream.str();
    file_name=file_name+"_"+rank_n+"_ref.msh";

    return file_name;
}

void mesh_parallel::fill_buffer_elements(const int& index){
  for(vector<int>::iterator j=i_buffer.begin();j!=i_buffer.end();j+=6){
    *j=my_elements.v((j-i_buffer.begin())/6+static_cast<int>(elmdist[index]))[0];
    *(j+1)=my_elements.v((j-i_buffer.begin())/6+static_cast<int>(elmdist[index]))[1];
    *(j+2)=my_elements.v((j-i_buffer.begin())/6+static_cast<int>(elmdist[index]))[2];
    *(j+3)=my_elements.v((j-i_buffer.begin())/6+static_cast<int>(elmdist[index]))[3];
    *(j+4)=my_elements.i((j-i_buffer.begin())/6+static_cast<int>(elmdist[index]));
    *(j+5)=index;
  }
}

void mesh_parallel::unpack_buffer_elements(const int& pos){
  for(vector<int>::iterator j=i_buffer.begin();j!=i_buffer.end();j+=6){
      my_elements.set_v(*j,0,pos+(j-i_buffer.begin())/6);
      my_elements.set_v(*(j+1),1,pos+(j-i_buffer.begin())/6);
      my_elements.set_v(*(j+2),2,pos+(j-i_buffer.begin())/6);
      my_elements.set_v(*(j+3),3,pos+(j-i_buffer.begin())/6);
      my_elements.set_i(*(j+4),pos+(j-i_buffer.begin())/6);
      my_elements.set_tag(*(j+5),pos+(j-i_buffer.begin())/6);
  }
}

void mesh_parallel::fill_buffer_ghosts(const int& index){
  for(vector<int>::iterator j=i_buffer.begin();j!=i_buffer.end();j+=6){
    *j=my_elements.v(ghosts[index][(j-i_buffer.begin())/6])[0];
    *(j+1)=my_elements.v(ghosts[index][(j-i_buffer.begin())/6])[1];
    *(j+2)=my_elements.v(ghosts[index][(j-i_buffer.begin())/6])[2];
    *(j+3)=my_elements.v(ghosts[index][(j-i_buffer.begin())/6])[3];
    *(j+4)=my_elements.i(ghosts[index][(j-i_buffer.begin())/6]);
    *(j+5)=find_pos_elmdist(ghosts[index][(j-i_buffer.begin())/6]);
  }
}

void mesh_parallel::fill_buffer_nodes(void){
  for(vector<int>::iterator j=i_buffer.begin();j!=i_buffer.end();++j){
    *j=my_nodes.i(j-i_buffer.begin());
    d_buffer[(j-i_buffer.begin())*2]=my_nodes.x(j-i_buffer.begin());
    d_buffer[(j-i_buffer.begin())*2+1]=my_nodes.y(j-i_buffer.begin());
  }
}

void mesh_parallel::unpack_buffer_nodes(void){
  for(vector<int>::iterator j=i_buffer.begin();j!=i_buffer.end();++j){
      my_nodes.set_x(d_buffer[(j-i_buffer.begin())*2],j-i_buffer.begin());
      my_nodes.set_y(d_buffer[(j-i_buffer.begin())*2+1],j-i_buffer.begin());
      my_nodes.set_i(*j,j-i_buffer.begin());
    }
}

void mesh_parallel::fill_buffer_lines(void){
  for(vector<int>::iterator j=i_buffer.begin();j!=i_buffer.end();j+=3){
    *j=my_lines.p((j-i_buffer.begin())/3)[0];
    *(j+1)=my_lines.p((j-i_buffer.begin())/3)[1];
    *(j+2)=my_lines.i((j-i_buffer.begin())/3);
  }
}

void mesh_parallel::unpack_buffer_lines(void){
  for(vector<int>::iterator j=i_buffer.begin();j!=i_buffer.end();j+=3){
      my_lines.set_p(*j,*(j+1),(j-i_buffer.begin())/3);
      my_lines.set_i(*(j+2),(j-i_buffer.begin())/3);
  }
}

void mesh_parallel::search_ghosts(void){

  ghosts.resize(static_cast<unsigned>(num_CPUs));

  for(int j=0;j!=num_CPUs;++j){
    for(int i=elmdist[j];i!=elmdist[j+1];++i){
      for(int k=0;k!=my_elements.num_neighbors(i);++k){
        if(((my_elements.n(i))[k]>-1)&&(my_elements.n(i)[k]<elmdist[j]||(my_elements.n(i))[k]>(elmdist[j+1]-1))){
          if(find(ghosts[j].begin(),ghosts[j].end(),(my_elements.n(i))[k])==ghosts[j].end())
            ghosts[j].push_back((my_elements.n(i))[k]);
        }
      }
    }
  }
}

void mesh_parallel::set_ghosts_elm_proc(void){
  ghosts.resize(static_cast<unsigned>(num_CPUs));
  elm_proc.resize(static_cast<unsigned>(num_CPUs));

  for(int j=0;j!=num_CPUs;++j){
    for(vector<long int>::iterator i=partition.begin();i!=partition.end();++i){
      if(j==*i) elm_proc[j].push_back(i-partition.begin());
    }
  }
  for(int j=0;j!=num_CPUs;++j){
    for(vector<int>::iterator i=elm_proc[j].begin();i!=elm_proc[j].end();++i){
      for(int k=0;k!=my_elements.num_neighbors(*i);++k){
        if((my_elements.n(*i)[k]>-1)){
          if(find(elm_proc[j].begin(),elm_proc[j].end(),my_elements.n(*i)[k])==elm_proc[j].end()){
            if(find(ghosts[j].begin(),ghosts[j].end(),(my_elements.n(*i))[k])==ghosts[j].end())
              ghosts[j].push_back((my_elements.n(*i))[k]);
          }
        }
      }
    }
  }
}

void mesh_parallel::pre_partitioning(void){

  // nodes broadcasting
  temp_dim=my_nodes.num_nodes();
  MPI_Bcast(&temp_dim,1,MPI_INT,0,MPI_COMM_WORLD);
  if(my_rank!=0) my_nodes.resize(temp_dim);
  d_buffer.resize(static_cast<unsigned>(2*my_nodes.num_nodes()));
  i_buffer.resize(static_cast<unsigned>(my_nodes.num_nodes()));
  if(my_rank==0) fill_buffer_nodes();
  MPI_Bcast(&(d_buffer[0]),2*my_nodes.num_nodes(),MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&(i_buffer[0]),my_nodes.num_nodes(),MPI_INT,0,MPI_COMM_WORLD);
  if(my_rank!=0) unpack_buffer_nodes();

  // elmdist creation and broadcasting
  elmdist.resize(static_cast<unsigned>(num_CPUs+1));
  if(my_rank==0){
    int num_elements_on_each_CPU(my_elements.num_elements()/num_CPUs),num_exciding_elements(my_elements.num_elements()%num_CPUs),temp(0);
    elmdist[0]=0;
    for(vector<long int>::iterator i=(elmdist.begin()+1);i!=elmdist.end();++i){
      temp=temp+num_elements_on_each_CPU;
      if(static_cast<int>((i-elmdist.begin()))<(num_exciding_elements+1)) ++temp;
      *i=temp;
    }
  }
  MPI_Bcast(&(elmdist[0]),num_CPUs+1,MPI_LONG,0,MPI_COMM_WORLD);

  // elements scattering for rank>0
  if(my_rank==0){
    for(int j=1;j!=num_CPUs;++j){
      i_buffer.resize(static_cast<unsigned>(6*(elmdist[j+1]-elmdist[j])));
      fill_buffer_elements(j);
      MPI_Send(&(i_buffer[0]),static_cast<int>(6*(elmdist[j+1]-elmdist[j])),MPI_INT,j,1,MPI_COMM_WORLD);
    }
  }
  else{
    my_elements.resize(static_cast<int>(elmdist[my_rank+1]-elmdist[my_rank]));
    i_buffer.resize(static_cast<unsigned>(6*my_elements.num_elements()));
    MPI_Recv(&(i_buffer[0]),6*my_elements.num_elements(),MPI_INT,0,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    unpack_buffer_elements(0);
    }

  // ghost scattering dimension
  if(my_rank==0){
    search_ghosts();
    for(vector<vector<int> >::iterator j=(ghosts.begin()+1);j!=ghosts.end();++j){
      temp_dim=static_cast<int>(j->size());
      MPI_Send(&temp_dim,1,MPI_INT,(j-ghosts.begin()),2,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Recv(&temp_dim,1,MPI_INT,0,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    my_elements.resize(my_elements.num_elements()+temp_dim);
    my_elements.set_num_ghosts(temp_dim);
  }

  // ghost elements scattering
  if(my_rank==0){
    for(int j=1;j!=num_CPUs;++j){
        i_buffer.resize(6*ghosts[j].size());
        fill_buffer_ghosts(j);
        MPI_Send(&(i_buffer[0]),static_cast<int>(6*ghosts[j].size()),MPI_INT,j,3,MPI_COMM_WORLD);
    }
  }
  else{
    i_buffer.resize(6*static_cast<unsigned>(my_elements.num_ghosts()));
    MPI_Recv(&(i_buffer[0]),6*my_elements.num_ghosts(),MPI_INT,0,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    unpack_buffer_elements(num_elements());
  }

  // lines broadcast
  if(my_rank==0) temp_dim=my_lines.num_lines();
  MPI_Bcast(&temp_dim,1,MPI_INT,0,MPI_COMM_WORLD);
  if(my_rank!=0) my_lines.resize(temp_dim);
  i_buffer.resize(3*static_cast<unsigned>(my_lines.num_lines()));
  if(my_rank==0) fill_buffer_lines();
  MPI_Bcast(&(i_buffer[0]),3*my_lines.num_lines(),MPI_INT,0,MPI_COMM_WORLD);
  if(my_rank!=0) unpack_buffer_lines();

  // resize my_elements on processor 0
  if(my_rank==0){
    sort(ghosts[0].begin(),ghosts[0].end());
    for(int i=0;i!=static_cast<int>(ghosts[0].size());++i){
      my_elements.copy(ghosts[0][i],elmdist[1]+i);
      my_elements.set_tag(find_pos_elmdist(ghosts[0][i]),elmdist[1]+i);
    }
    my_elements.resize(elmdist[1]+static_cast<int>(ghosts[0].size()));
    my_elements.set_num_ghosts(static_cast<int>(ghosts[0].size()));
    for(int j=0;j!=num_elements();++j) my_elements.set_tag(0,j);
    my_elements.free_neighbors();
  }
  free_unused_pre();
}

void mesh_parallel::create_ghosts_mapping(void){

  my_ghosts_comm->set_size(2);
  vector<double> values(2*num_elements(),0);
  for(int j=0;j!=(2*num_elements());j+=2){
    values[j]=my_rank;
    values[j+1]=j/2;
  }
  my_ghosts_comm->set_global_vec(values);

  ghosts_mapping.resize(static_cast<unsigned>(2*num_ghosts()));
  for(int j=0;j!=2*num_ghosts();++j){
    ghosts_mapping[j]=static_cast<int>(my_ghosts_comm->get_ghost_value(j));
    my_elements.set_tag(get_ghosts_rank(j/2),j/2);
  }
}

void mesh_parallel::set_options(const int& ncon_,const int& ncommonnodes_,const float& ubvecc_){
  ncon=ncon_;
  ncommonnodes=ncommonnodes_;
  ubvecc=ubvecc_;
}

void mesh_parallel::create_par_send_buffers(void){
  int dim_temp_ghosts(0),cont(0);
  send_elements_count.resize(static_cast<unsigned>(num_CPUs));
  send_ghosts_count.resize(static_cast<unsigned>(num_CPUs));
  for(int j=0;j!=num_CPUs;++j){
    send_elements_count[j]=static_cast<int>(5*elm_proc[j].size());
    send_ghosts_count[j]=static_cast<int>(5*ghosts[j].size());
    dim_temp_ghosts+=send_ghosts_count[j];
  }
  send_elements_buffer.resize(5*static_cast<unsigned>(my_elements.num_elements()));
  for(vector<vector<int> >::iterator j=elm_proc.begin();j!=elm_proc.end();++j){
    for(vector<int>::iterator k=j->begin();k!=j->end();++k){
      send_elements_buffer[5*cont]=my_elements.i(*k);
      send_elements_buffer[5*cont+1]=my_elements.v(*k)[0];
      send_elements_buffer[5*cont+2]=my_elements.v(*k)[1];
      send_elements_buffer[5*cont+3]=my_elements.v(*k)[2];
      send_elements_buffer[5*cont+4]=my_elements.v(*k)[3];
      ++cont;
    }
  }
  send_ghosts_buffer.resize(static_cast<unsigned>(dim_temp_ghosts));
  cont=0;
  for(vector<vector<int> >::iterator j=ghosts.begin();j!=ghosts.end();++j){
      for(vector<int>::iterator k=j->begin();k!=j->end();++k){
        send_ghosts_buffer[5*cont]=my_elements.i(*k);
        send_ghosts_buffer[5*cont+1]=my_elements.v(*k)[0];
        send_ghosts_buffer[5*cont+2]=my_elements.v(*k)[1];
        send_ghosts_buffer[5*cont+3]=my_elements.v(*k)[2];
        send_ghosts_buffer[5*cont+4]=my_elements.v(*k)[3];
        ++cont;
      }
    }
}

void mesh_parallel::unpack_par_elements_buffer(const int& num_new_elements){
  int old_size(my_elements.num_elements());
  my_elements.resize(old_size+num_new_elements);
  for(vector<int>::iterator j=receive_elements_buffer.begin();j!=receive_elements_buffer.end();++j){
    if((j-receive_elements_buffer.begin())%5==0){
      my_elements.set_i(*j,old_size+(j-receive_elements_buffer.begin())/5);
      my_elements.set_tag(my_rank,old_size+(j-receive_elements_buffer.begin())/5);
    }
    if((j-receive_elements_buffer.begin())%5==1) my_elements.set_v(*j,0,old_size+(j-receive_elements_buffer.begin())/5);
    if((j-receive_elements_buffer.begin())%5==2) my_elements.set_v(*j,1,old_size+(j-receive_elements_buffer.begin())/5);
    if((j-receive_elements_buffer.begin())%5==3) my_elements.set_v(*j,2,old_size+(j-receive_elements_buffer.begin())/5);
    if((j-receive_elements_buffer.begin())%5==4) my_elements.set_v(*j,3,old_size+(j-receive_elements_buffer.begin())/5);
  }
}

void mesh_parallel::unpack_par_ghosts_buffer(void){
  element new_ghost;
  bool flag=false;
  for(vector<int>::iterator k=receive_ghosts_buffer.begin();k!=receive_ghosts_buffer.end();k+=5){
    for(int j=0;(j!=(my_elements.num_elements()+my_elements.num_ghosts()))&&(flag==false);++j){
      if(my_elements.i(j)==*k) flag=true;
    }
    if(flag==false){
      new_ghost.set_i(*k);
      new_ghost.set_v(*(k+1),0);
      new_ghost.set_v(*(k+2),1);
      new_ghost.set_v(*(k+3),2);
      new_ghost.set_v(*(k+4),3);
      new_ghost.set_tag(0);
      my_elements.add_ghost(new_ghost);
    }
    flag=false;
  }
}

void mesh_parallel::partitioning(){

  ubvec.resize(static_cast<unsigned>(ncon));
  tpwgts.resize(static_cast<unsigned>(ncon*num_CPUs));

  for(vector<float>::iterator i=ubvec.begin();i!=ubvec.end();++i) *i=ubvecc;
  for(vector<float>::iterator i=tpwgts.begin();i!=tpwgts.end();++i) *i=1/(static_cast<float>(num_CPUs));

  partition.resize(static_cast<unsigned>(my_elements.num_elements()),my_rank);

  input_parmetis_conversion();

  MPI_Comm *comm_ptr=new MPI_Comm(MPI_COMM_WORLD);
  long int num_CPUs_long=static_cast<long int>(num_CPUs);
  long int temp(2),temp1(0);

  ParMETIS_V3_PartMeshKway(&(elmdist[0]),&(eptr[0]),&(eind[0]),&(weights[0]),&temp,&temp1,&ncon,&ncommonnodes,
    &(num_CPUs_long),&(tpwgts[0]),&(ubvec[0]),&(options[0]),&edgecut,&(partition[0]),comm_ptr);

  partitioning_move_elements();
}

void mesh_parallel::partitioning_move_elements(void){
  set_ghosts_elm_proc();
  create_par_send_buffers();
  my_elements.free();
  int temp_size(0),cont(0);
  vector<int> displs(num_CPUs);
  for(vector<int>::iterator i=displs.begin();i!=displs.end();++i){
    *i=cont;
    cont+=send_elements_count[i-displs.begin()];
  }
  for(int j=0;j!=num_CPUs;++j){
    MPI_Scatter(&(send_elements_count[0]),1,MPI_INT,&temp_size,1,MPI_INT,j,MPI_COMM_WORLD);
    receive_elements_buffer.resize(static_cast<unsigned>(temp_size));
    MPI_Scatterv(&(send_elements_buffer[0]),&(send_elements_count[0]),&(displs[0]),MPI_INT,&(receive_elements_buffer[0]),temp_size,MPI_INT,j,MPI_COMM_WORLD);
    unpack_par_elements_buffer(temp_size/5);
  }
  cont=0;
  for(vector<int>::iterator i=displs.begin();i!=displs.end();++i){
      *i=cont;
      cont+=send_ghosts_count[i-displs.begin()];
  }
  for(int j=0;j!=num_CPUs;++j){
    MPI_Scatter(&(send_ghosts_count[0]),1,MPI_INT,&temp_size,1,MPI_INT,j,MPI_COMM_WORLD);
    receive_ghosts_buffer.resize(static_cast<unsigned>(temp_size));
    MPI_Scatterv(&(send_ghosts_buffer[0]),&(send_ghosts_count[0]),&(displs[0]),MPI_INT,&(receive_ghosts_buffer[0]),temp_size,MPI_INT,j,MPI_COMM_WORLD);
    unpack_par_ghosts_buffer();
  }
  free_unused_par();
}

void mesh_parallel::compute_elmdist(void){
  elmdist.resize(static_cast<unsigned>(num_CPUs+1),0);
  long int temp(my_elements.num_elements());
  MPI_Gather(&temp,1,MPI_LONG,&elmdist[1],1,MPI_LONG,0,MPI_COMM_WORLD);
  if(my_rank==0){
    for(vector<long int>::iterator iter=(elmdist.begin()+2);iter!=elmdist.end();++iter) *iter+=*(iter-1);
  }
  MPI_Bcast(&(elmdist[0]),num_CPUs+1,MPI_LONG,0,MPI_COMM_WORLD);
}

void mesh_parallel::refinement(void){

  long int temp1(0), temp(2);
  MPI_Comm *comm_ptr=new MPI_Comm(MPI_COMM_WORLD);
  long int num_CPUs_long(num_CPUs);

  partition.resize(static_cast<unsigned>(my_elements.num_elements()),my_rank);
  input_parmetis_conversion();
  ParMETIS_V3_Mesh2Dual(&(elmdist[0]),&(eptr[0]),&(eind[0]),&temp1,&ncommonnodes,&xadj_ptr,&adjncy_ptr,comm_ptr);
  ParMETIS_V3_RefineKway(&(elmdist[0]),xadj_ptr,adjncy_ptr,&(weights[0]),NULL,&temp,&temp1,&ncon,&num_CPUs_long,&(tpwgts[0]),&(ubvec[0]),&(options[0]),&edgecut,&(partition[0]),comm_ptr);
  partitioning_move_elements();
  eval_topology();
  compute_elmdist();
}

void mesh_parallel::input_parmetis_conversion(void){

  eptr.resize(static_cast<unsigned>(my_elements.num_elements())+1);
  eind.resize(static_cast<unsigned>(my_elements.num_elements())*4);

  for(vector<long int>::iterator i=eptr.begin();i!=eptr.end();++i){
    int pos(i-eptr.begin());
    *i=pos*4;
    if(pos<(my_elements.num_elements())){
      for(int j=0;j!=4;++j) eind[pos*4+j]=static_cast<long int>(my_elements.v(pos)[j]);
    }
  }
}

void mesh_parallel::set_periodicity(const int& starting_line,const int& ending_line,const int& cor_starting_line,const int& cor_ending_line){

  set_periodicity_mesh(starting_line,ending_line,cor_starting_line,cor_ending_line);

  my_ghosts_comm->set_size(4);
  vector<double> nodes(num_elements()*4);
  for(int j=0;j!=num_elements();++j){
    nodes[j*4]=my_elements.v(j)[0];
    nodes[j*4+1]=my_elements.v(j)[1];
    nodes[j*4+2]=my_elements.v(j)[2];
    nodes[j*4+3]=my_elements.v(j)[3];
  }
  my_ghosts_comm->set_global_vec(nodes);

  for(int j=num_elements();j!=elements_size();++j){
    my_elements.set_v(static_cast<int>(my_ghosts_comm->get_ghost_value((j-num_elements())*4)),0,j);
    my_elements.set_v(static_cast<int>(my_ghosts_comm->get_ghost_value((j-num_elements())*4+1)),1,j);
    my_elements.set_v(static_cast<int>(my_ghosts_comm->get_ghost_value((j-num_elements())*4+2)),2,j);
    my_elements.set_v(static_cast<int>(my_ghosts_comm->get_ghost_value((j-num_elements())*4+3)),3,j);
  }

}
