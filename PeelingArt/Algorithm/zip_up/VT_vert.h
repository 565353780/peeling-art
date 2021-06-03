#pragma once

class VT_vert
{
public:
	VT_vert()
	{
		deepth_ = -1;
		index_global_ = -1;
		index_local_ = -1;
		father_ = 0;
		degree_ = 0;
		sons_.resize(0);
	}

	/*从树上脱离出去*/
	void disconnect()
	{
		deepth_ = 0;
		father_ = 0;
		degree_ = 0;
		sons_.resize(0);
		
	}

	~VT_vert(){}

	int getDeepth()
	{
		if (deepth_ != -1)
		{
			return deepth_;
		}
		else
		{
			degree_++;
			(father_->degree_)++;
			return deepth_ = father_->getDeepth() + 1;
		}
	}

	int deepth_;

	int index_global_;

	int index_local_;

	int degree_;

	VT_vert *father_;

	std::vector<VT_vert *> sons_;
};
