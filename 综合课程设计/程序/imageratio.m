function cr = imageratio(f1, f2)
%��������ͼ��ѹ����
% error(nargchk(2, 2, nargin));
cr = bytes(f1)/bytes(f2);
disp('ѹ��ǰͼ������������������:');
bytes(f1)
disp('ѹ��ǰͼ������������������:');
bytes(f2)
