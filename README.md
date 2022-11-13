# PEN_PPsystem
PEN toolboxに加えてPredator Prey System(PP)も使用できる化学反応ネットワーク設計・シミュレーションシステムを実装した。
先行研究としてPENのみのシミュレーションやPEN＋PPの常微分方程式構築があったが、PEN+PPのシミュレーションツールは今まで実装されていなかった。
この開発システムのメリットとしてはユーザーが化学反応ネットワークを入力しただけで簡単に、反応タイプに応じた反応中間体の生成と、微分方程式に基づいたシミュレーションを行うことができる。すべてのコードはPythonで作成した。

今回OligaterとBistable　Switchのモデルを化学反応ネットワークとして入力した。
結果、入力の際設定した反応タイプ(PEN,PP)に応じて反応中間体を自動生成し、全ての反応の常微分方程式を解きシミュレーションまで一気通貫で行われた。実行結果は以下の通りだ。

・Oligater



・Bistable Switch(双安定性スイッチ）->PEN+PP
<img width="523" alt="Bistable" src="https://user-images.githubusercontent.com/93179388/201527889-b656225d-18f0-4e7a-963c-d90d3c071be7.png">
