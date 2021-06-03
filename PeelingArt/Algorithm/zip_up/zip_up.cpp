#include "zip_up.h"

#include <fstream>

using namespace std;

#define RG(K, A, B) for(int (K) = (A); (K) < (B); (K)++)

ZipUp::ZipUp()
{

}

ZipUp::~ZipUp()
{
}

void ZipUp::getBoundaryIndex()
{
	boundary_vertex_index_.resize(0);

	auto heit = mesh_.halfedges_begin();
	while (!mesh_.is_boundary(*heit))
		heit++;
	auto he_start = *heit;
	auto he_it = he_start;
	do
	{
		he_it = mesh_.next_halfedge_handle(he_it);
		boundary_vertex_index_.push_back(mesh_.to_vertex_handle(he_it).idx());
	} while (he_it != he_start);
}

void ZipUp::getBoundaryAngle()
{
	/*run after getBoundaryIndex*/
	int bvert_num = boundary_vertex_index_.size();

	boundary_angle_.resize(bvert_num);

	for (int k(0); k < bvert_num; k++)
	{
		int kp = (k + bvert_num - 1) % bvert_num;
		int kn = (k + 1) % bvert_num;

		OpenMesh::Vec3d xp  = vertex_list[boundary_vertex_index_[kp]];
		OpenMesh::Vec3d x   = vertex_list[boundary_vertex_index_[k]];
		OpenMesh::Vec3d xn  = vertex_list[boundary_vertex_index_[kn]];

		xp -= x;
		xn -= x;

		float cosin = OpenMesh::dot(xp, xn) / (xp.length()* xn.length());
		float ang = acos(cosin);
		boundary_angle_[k] = ang;
	}
}

void ZipUp::Prim( bool plan_B, bool is_uniform )
{
	/*run after getBoundaryIndex*/

	/*局部类*/
	class DataUpdater
	{
	public:

		~DataUpdater(){}

		DataUpdater(vector<float> &ml, vector<float> &ad, vector<float> &df, 
			vector<bool> &isns, vector<OpenMesh::Vec3d> &vl, vector<int> &bi, 
			bool plan_B, bool is_uniform)
			:ml_(ml),ad_(ad),df_(df),isns_(isns),vl_(vl),
			bi_(bi),plan_B_(plan_B),is_uniform_(is_uniform){}

		void updateAt(int kp, int k, int kn)
		{
			OpenMesh::Vec3d xp	= vl_[bi_[kp]];
			OpenMesh::Vec3d x = vl_[bi_[k]];
			OpenMesh::Vec3d xn = vl_[bi_[kn]];

			OpenMesh::Vec3d xpn = xn - xp;
			xp -= x;
			xn -= x;

			isns_[k] = (xp.length() > xn.length());

			if (! plan_B_)
			{
				/*plan_A*/
				ml_[k] = abs(xp.length() - xn.length());
			}
			else
			{
				/*plan_B*/
				ml_[k] = max(xp.length(), xn.length());
			}

			ad_[k] = xpn.length();
			df_[k] = ml_[k] - ad_[k];

			if (is_uniform_)
			{
				df_[k] /= xp.length() + xn.length() + xpn.length();
			}
		}

		/*第k个边界点关联边界边的长度最长值*/
		vector<float> &ml_;

		/*第k个边界点相邻边界点之间的距离*/
		vector<float> &ad_;

		/*上边两个的差*/
		vector<float> &df_;

		/*第k个边界点关联边界边的长度最长值*/
		vector<bool> &isns_;

		vector<OpenMesh::Vec3d> &vl_;

		vector<int> &bi_;

		/*是否采用plan_B解决方案*/
		bool plan_B_;

		/*是否归一化*/
		bool is_uniform_;
	};

	/*第k个边界点关联边界边的长度最长值*/
	vector<float> inc_max_length;

	/*第k个边界点相邻边界点之间的距离*/
	vector<float> adj_distance;

	/*上边两个的差*/
	vector<float> diff;

	/*第k个边界点关联边界边的长度最长值*/
	vector<bool> is_next_smaller;

	int bvert_num = boundary_vertex_index_.size();

	inc_max_length.resize(bvert_num);
	is_next_smaller.resize(bvert_num);
	adj_distance.resize(bvert_num);
	diff.resize(bvert_num);

	/*树结构*/
	tree_edges_.resize(2, bvert_num - 1);

	tvert_list_.resize(bvert_num);
	RG(k,0,bvert_num)
	{
		tvert_list_[k].index_global_ = boundary_vertex_index_[k];
		tvert_list_[k].index_local_ = k;
	}
	DataUpdater DU(inc_max_length, adj_distance, diff, is_next_smaller, 
		vertex_list, boundary_vertex_index_, plan_B, is_uniform);

	for (int k(0); k < bvert_num; k++)
	{
		int kp = (k + bvert_num - 1) % bvert_num;
		int kn = (k + 1) % bvert_num;

		DU.updateAt(kp, k, kn);
	}

	vector<int> ind;
	ind.resize(bvert_num);

	RG(i,0,bvert_num)
	{
		ind[i] = i;
	}

	RG(i,0,bvert_num - 2)
	{
		int inum = ind.size();

		float max_diff(-100.0f);
		int k0(0);
		RG(k,0,inum)
		{
			if (diff[ind[k]] > max_diff)
			{
				max_diff = diff[ind[k]];
				k0 = k;
			}
		}

		int kn	= (k0 + 1) % inum;
		int knn	= (k0 + 2) % inum;
		int kp	= (k0 + inum - 1) % inum;
		int kpp	= (k0 + inum - 2) % inum;


		if (is_next_smaller[ind[k0]])
		{
			tree_edges_(0, i) = ind[k0];
			tree_edges_(1, i) = ind[kn];
			tvert_list_[ind[k0]].father_ = &(tvert_list_[ind[kn]]);
			tvert_list_[ind[kn]].sons_.push_back(&(tvert_list_[ind[k0]]));
		}
		else
		{
			tree_edges_(0, i) = ind[k0];
			tree_edges_(1, i) = ind[kp];
			tvert_list_[ind[k0]].father_ = &(tvert_list_[ind[kp]]);
			tvert_list_[ind[kp]].sons_.push_back(&(tvert_list_[ind[k0]]));
		}

		/*更新这几个向量*/
		DU.updateAt(ind[kp], ind[kn], ind[knn]);
		DU.updateAt(ind[kpp], ind[kp], ind[kn]);

		ind.erase(ind.begin() + k0);
 	}

	tree_edges_(0, bvert_num - 2) = ind[0];
	tree_edges_(1, bvert_num - 2) = ind[1];

	tvert_list_[ind[1]].father_ = &(tvert_list_[ind[0]]);
	tvert_list_[ind[0]].deepth_ = 0;

	tree_root_ = &(tvert_list_[ind[0]]);
	tvert_list_[ind[0]].father_ = NULL;
	tvert_list_[ind[0]].sons_.push_back(&(tvert_list_[ind[1]]));

	/*生成树结构*/
	RG(k,0,bvert_num)
		tvert_list_[k].getDeepth();

	bedge_structure_.resize(bvert_num);

	RG(k,0,bvert_num)
	{
		VT_vert *v0 = & tvert_list_[k];
		VT_vert *vn = & tvert_list_[(k + 1) % bvert_num];

		vector<int> &ks = bedge_structure_[k];

		auto itr = ks.begin();
		while(v0->index_global_ != vn->index_global_)
		{
			if (v0->deepth_ >= vn->deepth_)
			{
				v0 = v0->father_;
				if (v0->index_global_ != vn->index_global_)
				{
					itr = ks.insert(itr, v0->index_global_);
					itr++;
				}
				
			}
			else if(v0->deepth_ < vn->deepth_)
			{
				vn = vn->father_;
				if (v0->index_global_ != vn->index_global_)
				{
					itr = ks.insert(itr, vn->index_global_);
				}
			}
		}
	}


}

void ZipUp::writeBoundaryIndex( string a )
{
	ofstream file(a, ios::out);

	int bvert_num = boundary_vertex_index_.size();

	RG(k,0,bvert_num)
	{
		file << boundary_vertex_index_[k] << endl;
	}

	file.close();
}

void ZipUp::writeBoundaryStructure( string a )
{
	/*run after Prim*/

	ofstream file(a, ios::out);

	int bvert_num = boundary_vertex_index_.size();

	RG(k,0,bvert_num)
	{
		int kn = (k + 1) % bvert_num;
		file << boundary_vertex_index_[k] << ' ';

		vector<int> &ks = bedge_structure_[k];

		RG(i,0,ks.size())
		{
			file << ks[i] << ' ';
		}

		file << boundary_vertex_index_[kn] << endl;
	}

	file.close();
}

void ZipUp::writeTreeEdge( std::string a )
{
	/*run after Prim*/

	ofstream file(a, ios::out);

	int tedge_num = tree_edges_.cols();

	RG(k,0,tedge_num - 1)
	{
		file << boundary_vertex_index_[tree_edges_(0, k)] << ' ';
	}

	file << boundary_vertex_index_[tree_edges_(0, tedge_num - 1)] << endl;

	RG(k,0,tedge_num - 1)
	{
		file << boundary_vertex_index_[tree_edges_(1, k)] << ' ';
	}

	file << boundary_vertex_index_[tree_edges_(1, tedge_num - 1)];

	file.close();
}

void ZipUp::evaluateComponent()
{
	int bvert_num = tvert_list_.size();

	tree_component_.resize(0);
	component_length_.resize(0);

	vector<int> path_temp;
	float length_temp;

	/*修正根节点*/
	if(tree_root_->degree_ == 2)
	{
		do 
		{
			VT_vert *eldest_son = tree_root_->sons_[0];
			tree_root_->father_ = eldest_son;
			tree_root_->sons_.erase(tree_root_->sons_.begin());
			eldest_son->father_ = NULL;
			eldest_son->sons_.push_back(tree_root_);

			/*传位于嫡长子*/
			tree_root_ = eldest_son;
		} while (tree_root_->degree_ == 2);
	}

	RG(k,0,bvert_num)
	{
		VT_vert *vert_temp = &tvert_list_[k];
		if (vert_temp->degree_ != 2 && vert_temp->degree_ != 0 && vert_temp->father_ != NULL)
		{
			path_temp.resize(0);
			length_temp = 0;

			path_temp.push_back(vert_temp->index_local_);

			while (vert_temp->father_ != NULL)
			{
				int ind1 = vert_temp->index_global_;
				vert_temp = vert_temp->father_;
				int ind2 = vert_temp->index_global_;

				path_temp.push_back(vert_temp->index_local_);
				length_temp += (vertex_list[ind1] - vertex_list[ind2]).length();

				if (vert_temp->degree_ != 2)
				{
					break;
				}
			}

			tree_component_.push_back(path_temp);
			component_length_.push_back(length_temp);
		}
		
	}
}

void ZipUp::cleanComponent( float tol )
{
	int k(0);

	while (true)
	{
		vector<float> c_length;
		vector<int> ind_list;

		c_length.resize(0);
		ind_list.resize(0);

		RG(i,0,tree_component_.size())
		{
			vector<int> &path_temp = tree_component_[i];
			if (tvert_list_[path_temp[0]].degree_ == 1)
			{
				c_length.push_back(component_length_[i]);
				ind_list.push_back(i);
			}
		}

		if (ind_list.empty())
		{
			break;
		}

		vector<float>::iterator mini = std::min_element(begin(c_length), end(c_length));
		k = distance(begin(c_length), mini);

		if (c_length[k] > tol)
		{
			break;
		}

		k = ind_list[k];

		vector<int> &comp_dlting = tree_component_[k];
		RG(i,0,comp_dlting.size() - 1)
		{
			tvert_list_[comp_dlting[i]].disconnect();
		}

		tvert_list_[comp_dlting.back()].degree_ --;
		vector<VT_vert *> &sons_temp = tvert_list_[comp_dlting.back()].sons_;
 		for(auto son_iter = sons_temp.begin();son_iter != sons_temp.end();son_iter ++)
		{
			if ((*son_iter)->index_local_ == comp_dlting[comp_dlting.size() - 2])
			{
				sons_temp.erase(son_iter);
				break;
			}
 		}
		evaluateComponent();
	}
}

void ZipUp::writeComponent( std::string a )
{
	ofstream file(a, ios::out);

	int cmp_num = tree_component_.size();

	RG(k,0,cmp_num)
	{
		vector<int> &ck = tree_component_[k];

		RG(i,0,ck.size() - 1)
		{
			file << tvert_list_[ck[i]].index_global_ << ' ';
		}

		file << tvert_list_[ck.back()].index_global_ << endl;
	}

	file.close();
}

void ZipUp::writeBoundaryPosition( std::string a )
{
	ofstream file(a, ios::out);

	int bvert_num = boundary_vertex_index_.size();

	RG(k,0,bvert_num)
	{
		OpenMesh::Vec3d &pt = vertex_list[boundary_vertex_index_[k]];

		file << pt[0] << ' ' << pt[1] << ' ' << pt[2] << endl;
	}

	file.close();
}
